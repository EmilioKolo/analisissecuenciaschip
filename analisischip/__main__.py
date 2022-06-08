
import os
from pyensembl import EnsemblRelease
from Bio import Entrez
import logging
from analisischip import seq_data

'''
Archivo donde ejecuto las pruebas
Corre con python -m analisischip

Tipos de error:
DEBUG -> 10
INFO -> 20
WARNING -> 30
ERROR -> 40
CRITICAL -> 50
(usar logging.debug, logging.info, logging.warning, logging.error, logging.critical)
* Solo mensajes de warning para arriba van a consola
'''

#################################### VARIABLES ####################################

Entrez.email = 'ekolomenski@gmail.com';
Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';
mm9 = EnsemblRelease(67, species='mouse');
hg19 = EnsemblRelease(102, species='human');
GRCm38 = EnsemblRelease(102, species='mouse');
dict_genomas = {'mm9':mm9, 'mouse':mm9, 'hg19':hg19, 'human':hg19, 'GRCm38':GRCm38, 'mouse102':GRCm38};

L_SU_NKX25 = ['GCAAGTG', 'GGAAGTG', 'GAAAGTG', 'ATAAGTG', 'GTAAGTG',
              'CTAAGTG', 'TCAAGTG', 'TGAAGTG', 'TAAAGTG', 'TTAAGTG'];
rango1500 = (-1000, 1000);
L_nombres_chip_mm9 = ['Dupays2015', 'He2011', 'vandenBoogard2012'];
L_nombres_chip_dupays = ['Dupays_GSM1087143_nkx_s1e14_peaks', 'Dupays_GSM1087144_nkx_s4e14_peaks'];
L_nombres_chip_hg19 = ['Anderson2018-GSE89457consensus'];

path_bed = '..\\0-Fuentes\\Papers ChIP-seq\\';

curr_path = os.path.dirname(__file__);

# Variables main()

genoma_main = 'mm9';
L_sitios_main = L_SU_NKX25;
rango_sitios_main = rango1500;
L_nombres_chip_main = L_nombres_chip_mm9;

#################################### FUNCIONES ####################################


def _main():
    # Funcion para probar funciones en ejecucion del archivo

    # Inicializo la clase
    NKX25_SU = seq_data('mm9', dict_genomas);
    ### PRUEBA
    NKX25_SU.agregar_secuencia('chr1', 2200, 2600);
    print(NKX25_SU.dict_range)
    print()
    NKX25_SU.agregar_secuencia('chr2', 2100, 2700);
    print(NKX25_SU.dict_range)
    print()
    NKX25_SU.agregar_secuencia('chr1', 2500, 2900);
    print(NKX25_SU.dict_range)
    print()
    NKX25_SU.agregar_secuencia('chr1', 3500, 3900);
    print(NKX25_SU.dict_range)
    print()
    NKX25_SU.agregar_secuencia('chr1', 1500, 1900);
    print(NKX25_SU.dict_range)
    print()
    #NKX25_SU.revisar_mult(carga_vacias=True);
    ###
    print(NKX25_SU.M_seq)

    return NKX25_SU

#################################### RESULTADOS ###################################

output_dump = [];

if __name__=='__main__':
    output_dump.append(_main());
