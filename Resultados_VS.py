#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import time
import pandas as pd
import matplotlib.pyplot as plt
from statistics import mean
import sys
import shutil


__author__ = "Juan Manuel Prieto"


def Grafico_Clusters_Pob(enegy_df, caso):
    Eje_X = enegy_df['mean energy'].tolist()
    Eje_Y = enegy_df['cluster'].tolist()

    fig, ax = plt.subplots()
    ax.scatter(Eje_X, Eje_Y, label='Docking')

    plt.grid()
    plt.title(f'{caso}')
    plt.legend(loc='upper right')
    plt.xlabel('Energy')
    plt.ylabel('% Population')
    plt.savefig(f'Cluster_Poblaciones_{caso}.jpg')
    shutil.copy(f'Cluster_Poblaciones_{caso}.jpg', f'{path_complete}/Resultados/')


def ranks(input_file):
    rank_value = 1
    ranks = {}

    for lines in input_file:
        if lines[0:6] == 'DOCKED':
            mark = (lines.split(' '))[1]
            if mark == 'USER':
                if 'Run' in lines:
                    rank = int(lines[22:25])
                    pdbqt_out = open('Run{}.pdbqt'.format(rank), 'w')
            if mark in ['REMARK', 'HETATM', 'ATOM', 'ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF', 'TER']:
                pdbqt_out.write(lines[8:])
            if mark in ['ENDMDL']:
                pdbqt_out.close()
        if 'RANKING' in lines:
            ranks[int(lines[13:20])] = rank_value
            rank_value = rank_value + 1
    pdbqt_out.close()

    # Renombrar el archivo
    for keys in ranks:
        os.rename('Run{}.pdbqt'.format(keys), 'Rank{}.pdbqt'.format(ranks[keys]))

    for j in range(1, rank_value):
        print(j)
        os.system('obabel -ipdbqt Rank{}.pdbqt -opdb -O rank{}.pdb'.format(j, j))


def Histograma(input_file):
    ### Histogramas ###
    enegy_df = pd.DataFrame(columns=['rank', 'lowest energy', 'run', 'mean energy', 'cluster'])
    for j in range(0, len(input_file)):
        if 'CLUSTERING HISTOGRAM' in input_file[j]:
            flag = True
            pos = j + 10
            while flag:
                if input_file[pos][0:5] != '_____':
                    data = input_file[pos].split('|')
                    valores = [int(data[0]), float(data[1]), int(data[2]), float(data[3]), int(data[4])]
                    temp_df = pd.DataFrame([valores], columns=enegy_df.columns)
                    enegy_df = pd.concat([enegy_df, temp_df], ignore_index=True)
                    pos = pos + 1
                else:
                    flag = False

    ## plot ##
    datos = enegy_df['cluster']
    etiquetas = enegy_df['rank'].astype(str)

    plt.bar(etiquetas, datos, color='skyblue', edgecolor='black', alpha=0.7)
    for i in range(len(datos)):
        plt.text(etiquetas[i], datos[i], str(datos[i]), ha='center', va='bottom')

    plt.xlabel('Ranks')
    plt.ylabel('% Population')
    plt.title(f'Histogram {caso}')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.show()
    plt.savefig(f'Histograma_{caso}.jpg')
    plt.close()
    shutil.copy(f'Histograma_{caso}.jpg', f'{path_complete}/Resultados/')
    return enegy_df


if __name__ == '__main__':
    Rec = '2FSI'
    path_complete = os.getcwd()
    Lista_Proteinas = os.listdir('Receptores')
    Ligandos = os.listdir('Ligands')

    try:
        os.mkdir('{}'.format('Resultados'))
    except FileExistsError:
        pass

    try:
        os.mkdir(f'{path_complete}/Resultados/DLG')
    except FileExistsError:
        pass

    try:
        os.mkdir(f'{path_complete}/Resultados/DLG_BIAS')
    except FileExistsError:
        pass

    stats_final = pd.DataFrame(columns=['Lig', 'Rank', 'Energy', 'Cluster', 'Case'])
    failed_cases = []  # Lista para los casos que fallaron

    for Prot in Lista_Proteinas:
        Prot = Prot.split('.')[0]
        for lig in Ligandos:
            lig = lig.split('.')[0]
            print(f'{Prot} {lig}')
            try:
                os.chdir('{}_{}'.format(Prot, lig))
                file_name = f'{lig}.dlg'
                caso = f"{Rec} {lig}"
                input_file = open('{}'.format(file_name), 'r').readlines()
                ranks(input_file)
                DF_Energia = Histograma(input_file)
                for k in range(0, DF_Energia.shape[0]):
                    stats_final.loc[len(stats_final.index)] = [lig, DF_Energia.iloc[k, 0], DF_Energia.iloc[k, 3], DF_Energia.iloc[k, 4], 'Normal']
                Grafico_Clusters_Pob(DF_Energia, caso)
                shutil.copy('rank1.pdb', f'{path_complete}/Resultados/{lig}_rank1.pdb')
                # Mover el archivo DLG a la carpeta DLG
                shutil.move(file_name, f'{path_complete}/Resultados/DLG/{file_name}')
                os.chdir(path_complete)
            except Exception as e:
                failed_cases.append(file_name)
                os.chdir(path_complete)
                print(f"Error processing {file_name}: {e}")

            try:
                lig_folder = lig + '_Bias'
                print(f'{Prot} {lig}')
                os.chdir('{}_{}'.format(Prot, lig_folder))
                file_name_bias = f'{lig}.dlg'
                caso = f"{Rec} {lig} Bias"
                input_file = open('{}'.format(file_name_bias), 'r').readlines()
                ranks(input_file)
                DF_Energia = Histograma(input_file)
                Grafico_Clusters_Pob(DF_Energia, caso)
                shutil.copy('rank1.pdb', f'{path_complete}/Resultados/{lig}_Bias_rank1.pdb')
                # Mover el archivo DLG_BIAS a la carpeta DLG_BIAS
                shutil.move(file_name_bias, f'{path_complete}/Resultados/DLG_BIAS/{file_name_bias}')
                os.chdir(path_complete)
            except Exception as e:
                failed_cases.append(file_name_bias)
                os.chdir(path_complete)
                print(f"Error processing {file_name_bias}: {e}")

    stats_final.to_csv(f'{path_complete}/Resultados/Energy_cluster_ranks.csv')

    # Guardar los casos que fallaron en un archivo de texto
    with open(f'{path_complete}/Resultados/failed_cases.txt', 'w') as f:
        for item in failed_cases:
            f.write("%s\n" % item)


