#!/usr/bin/env python

"""Alinear_PDB.py: Dado 2 PDB , alinea estrucutra "Bound" sobre "Apo"."""

__author__      = "Juan M. Prieto"
__copyright__   = "2023"
__version__ = "1.0"


import argparse
from Bio import PDB
from Bio.PDB import Superimposer

def calcular_alineamiento_rmsd(archivo_apo, archivo_bound, cadena_id='A'):
    # Cargar las estructuras desde los archivos PDB
    estructura1 = PDB.PDBParser().get_structure("estructura1", archivo_apo)
    estructura2 = PDB.PDBParser().get_structure("estructura2", archivo_bound)

    # Obtener el modelo y la cadena (ajustar según tus necesidades)
    modelo1 = estructura1[0][cadena_id]
    modelo2 = estructura2[0][cadena_id]

    # Obtener los residuos de la cadena A
    residuos_A = list(modelo1)

    # Obtener el número del primer y último aminoácido
    primer_aminoacido_A = residuos_A[0].id[1]
    ultimo_aminoacido_A = residuos_A[-1].id[1]

    # Obtener los residuos de la cadena B
    residuos_B = list(modelo2)

    # Obtener el número del primer y último aminoácido
    primer_aminoacido_B = residuos_B[0].id[1]
    ultimo_aminoacido_B = residuos_B[-1].id[1]

    # Calcular el inicio y fin del tramo
    inicio_tramo = max(primer_aminoacido_A, primer_aminoacido_B)
    fin_tramo = min(ultimo_aminoacido_A, ultimo_aminoacido_B)

    # Seleccionar los átomos del tramo para ambas estructuras
    seleccion1 = PDB.Selection.unfold_entities(modelo1, cadena_id)[inicio_tramo - 1:fin_tramo]
    seleccion2 = PDB.Selection.unfold_entities(modelo2, cadena_id)[inicio_tramo - 1:fin_tramo]

    # Crear el objeto Superimposer y calcular la transformación de alineamiento
    superimposer = Superimposer()
    superimposer.set_atoms(seleccion1, seleccion2)
    superimposer.apply(modelo2.get_atoms())

    # Calcular el RMSD
    rmsd = superimposer.rms

    # Guardar la estructura alineada en un nuevo archivo PDB
    io = PDB.PDBIO()
    io.set_structure(estructura2)
    output_pdb = f"{archivo_bound.split('.')[0]}_alig.pdb"
    io.save(output_pdb)

    return inicio_tramo, fin_tramo, rmsd

def main():
    parser = argparse.ArgumentParser(description='Realizar alineamiento estructural y calcular RMSD.')
    parser.add_argument('-A', '--Apo_Proteina', required=True, help='Archivo PDB de la proteína apo')
    parser.add_argument('-B', '--Bound_Proteina', required=True, help='Archivo PDB de la proteína bound')
    parser.add_argument('-C', '--Chain', default='A', help='Identificador de cadena (por defecto: A)')

    args = parser.parse_args()

    inicio_tramo, fin_tramo, rmsd = calcular_alineamiento_rmsd(args.Apo_Proteina, args.Bound_Proteina, args.Chain)

    # Guardar resultados en un archivo de texto
    output_file = f"{args.Bound_Proteina.split('.')[0]}_resultados.txt"
    with open(output_file, 'w') as f:
        f.write(f"Alineamiento entre: {args.Apo_Proteina.split('.')[0]} y {args.Bound_Proteina.split('.')[0]}\n")
        f.write(f"Inicio Alineamiento: {inicio_tramo}\n")
        f.write(f"Fin Alineamiento: {fin_tramo}\n")
        f.write(f"RMSD: {rmsd:.4f} Å\n")

    print(f"Resultados guardados en {output_file}")

if __name__ == "__main__":
    main()
