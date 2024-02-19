import requests
import json
import time
import os
from Bio import pairwise2
from Bio.Align import substitution_matrices


def check_pdb_for_ligand(pdb_id_to_check):
    
    peptide = ["ALA","LYS","VAL","ARG","GLY","CYS","TYR", "PHE","LEU","ILE","MET","SER","PRO","THR","HIS","ASN","GLN","GLU","ASP","TRP"]
    family_nad = ["NAD","NAI","NDC","TAP","ADJ","NAJ","NDO","ZID","CAN","NAP","NDP","CND","NAQ","NHD","DND","NAX","NHO","NAC","NDB","ODP","NAE","NBP","PAD","NAH","NDA","SND"]
    family_fad = ["FAD","FMN","6FA","FNS","FAA","MGD","FAB","RFL","FAE","FAS","FDA","FMA"]
    nucl =["AMP","GDP","UDP","ADP","GNP","UTP","ANP", "GTP","PSU","ATP","2GP","CMP","TMP","CDP","TDP","CTP","TTP","GMP","UMP"]
    mol_junk = ["UNK", "UNX" "BMA","FUC","MAN","POP","BOG","GAL","MES","PYR","C8E","GLC","MPD","SPM","CIT","GOL","MYR","TRS","CRY","HED","NAG","XYS","DTT","LDA","NGA","EPE","LI1","PEG","F6P","MAL","PG4"]

    Resultados = []

    url=f'https://www.ebi.ac.uk/pdbe/search/pdb/select?q=pdb_id:{pdb_id_to_check}&wt=json'


    try:
        # Realizar la solicitud a la API
        respuesta = requests.get(url, timeout=30)
        respuesta.raise_for_status()  # Lanza una excepción si la solicitud no fue exitosa

        # Analizar la respuesta en formato XML
        info_componente = json.loads(respuesta.text)
    except requests.exceptions.RequestException as e:
        print(f'Error al realizar la solicitud: {e}')
        return(Resultados , 0)
    
    
    Meta_Data = (dict(info_componente['response']))

    try:
        Num_chains = (Meta_Data['docs'][0]['number_of_protein_chains'])
    except KeyError:
        Num_chains = '-'

    try:
        type_assembly = (Meta_Data['docs'][0]['assembly_type'][0])
    except KeyError:
        type_assembly = '-'
    
    ## Tiene Ligando? 
    try:
        Ligands = (Meta_Data['docs'][0]['q_interacting_ligands'])
    except KeyError:
        Resultados.append([Num_chains,type_assembly,'-',0,0])
        Ligands = []


    for j in range(0,len(Ligands)):
        l = Ligands[j].split(':')[0].strip()
        if l in (peptide + nucl + mol_junk + family_fad + family_nad):
            continue  
        else:
            Peso, Num_Atomos = het_id(l)
            if Num_Atomos > 10:
                lig = []
                Resultados.append([Num_chains,type_assembly,l,Peso, Num_Atomos])

    Sec = (Meta_Data['docs'][0]['molecule_sequence'])
    return(Resultados , Sec)
    
def het_id(componente_id):
    info_componente = buscar_info_componente_quimico(componente_id)
    info_componente_atom = buscar_atomo_componente_quimico(componente_id)
    if (info_componente) and (info_componente_atom):
        Peso = (info_componente[componente_id][0]["weight"])
        Tamaño = (len(info_componente_atom[componente_id]))
        return(Peso , Tamaño)
    return(0,0)
    
def buscar_info_componente_quimico(component_id):
    # URL base de la API del RCSB PDB para el Chemical Components Dictionary (CCD)
    url_base = 'https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/'

    # Construir la URL para la búsqueda del componente químico por ID
    url = f'{url_base}{component_id}'

    if component_id in estructuras_inexistentes:
        return None
    try:
        # Realizar la solicitud a la API
        respuesta = requests.get(url, timeout=30)
        respuesta.raise_for_status()  # Lanza una excepción si la solicitud no fue exitosa

        # Analizar la respuesta en formato XML
        info_componente = json.loads(respuesta.text)

        return (info_componente)
    except requests.exceptions.RequestException as e:
        estructuras_inexistentes.append(component_id)
        print(f'Error al realizar la solicitud: {e}')
        return None


def buscar_atomo_componente_quimico(component_id):
    # URL base de la API del RCSB PDB para el Chemical Components Dictionary (CCD)
    url_base = 'https://www.ebi.ac.uk/pdbe/api/pdb/compound/atoms/'

    # Construir la URL para la búsqueda del componente químico por ID
    url = f'{url_base}{component_id}'

    if component_id in estructuras_inexistentes:
        return None

    try:
        # Realizar la solicitud a la API
        respuesta = requests.get(url, timeout=30)
        respuesta.raise_for_status()  # Lanza una excepción si la solicitud no fue exitosa

        # Analizar la respuesta en formato XML
        info_componente = json.loads(respuesta.text)

        return (info_componente)
    except requests.exceptions.RequestException as e:
        estructuras_inexistentes.append(component_id)
        print(f'Error al realizar la solicitud: {e}')
        return None

def Alinear_y_revisar(Sec_Apo, Sec_Bound):
    matrix = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.globalxx(Sec_Apo, Sec_Bound)
    Secuencia1 = (alignments[0].seqA)
    Secuencia2 = (alignments[0].seqB)
    Largo = len(Secuencia1)
    Dif = 0
    for j in range(0,len(Secuencia1)):
        if (Secuencia1[j]) != Secuencia2[j]:
            if (j > 10) and (j < (Largo - 10)):
                Dif = Dif + 1
    if Dif == 0:
        return(True)
    else:
        return(False)

if __name__ == '__main__':

    file_path = 'Clusters.txt'
    try:
        os.remove(file_path)
    except OSError as e:
        pass
    resultados = open("Clusters.txt", 'a')
    resultados.write('cluster,pdb,num_chain,assable,lig,weight,atoms\n')
    resultados.close()


    # Lee el archivo
    with open('clusters-by-entity-95-filtred-chains8.txt', 'r') as archivo:
        lineas = archivo.readlines()


    # Procesa cada línea
    entity_ids = [linea.strip().split() for linea in lineas]

    estructuras_inexistentes = []

    for j in range(0, len(entity_ids)):

        resultados = open("Clusters.txt", 'a')
        entity = entity_ids[j]

        PDB_con_ligandos = 0

        print(f'Cluster {j}')
        ### Busco s/lignado ###
        for pdb in entity:
            
            
            pdb_id_to_check = pdb[0:4].lower()  # Reemplazar con el ID del PDB que deseas verificar
            Datos_Ligandos , Secuencia = check_pdb_for_ligand(pdb_id_to_check)
            print(f'Cluster {j} Apo PDB: {pdb_id_to_check}')
            if len(Datos_Ligandos) == 0:
                entity.remove(pdb)
            elif (Datos_Ligandos[0][2] == '-') :
                resultados.write(f'{j},{pdb_id_to_check},{Datos_Ligandos[0][0]},{Datos_Ligandos[0][1]},{Datos_Ligandos[0][2]},{Datos_Ligandos[0][3]},{Datos_Ligandos[0][4]},0\n')
                Apo = pdb_id_to_check
                Sec_Apo = Secuencia
                entity.remove(pdb)
                break 
        ### Busco c/lignado ###
        for pdb in entity:
            pdb_id_to_check = pdb[0:4].lower()  # Reemplazar con el ID del PDB que deseas verificar
            Datos_Ligandos , Secuencia = check_pdb_for_ligand(pdb_id_to_check)
            
            print(f'Cluster {j} Bound PDB: {pdb_id_to_check} encontrados:{PDB_con_ligandos}')
            if len(Datos_Ligandos) == 0:
                entity.remove(pdb)
            elif (Datos_Ligandos[0][2] != '-') :
                ### agregar alineamiento ###
                Alineado = Alinear_y_revisar(Sec_Apo, Secuencia)
                if Alineado == True:
                    for k in range(0,len(Datos_Ligandos)):
                        resultados.write(f'{j},{pdb_id_to_check},{Datos_Ligandos[k][0]},{Datos_Ligandos[k][1]},{Datos_Ligandos[k][2]},{Datos_Ligandos[k][3]},{Datos_Ligandos[k][4]}\n')
                    PDB_con_ligandos = PDB_con_ligandos + 1
            if PDB_con_ligandos == 3:
                break
        resultados.close()
        time.sleep(30)


 


#Error al realizar la solicitud: HTTPSConnectionPool(host='www.ebi.ac.uk', port=443): Read timed out.