import logging
from analisischip import seq_data

'''
Archivo donde ejecuto las pruebas
Corre con python -m analisischip

Tipos de error:
INFO -> 10
DEBUG -> 20
WARNING -> 30
ERROR -> 40
CRITICAL -> 50
(usar logging.info, logging.debug, logging.warning, logging.error, logging.critical)
* Solo mensajes de warning para arriba van a consola
'''

# Dejar hasta terminar development
logging.basicConfig(level=logging.DEBUG);

output_dump = [];

if __name__=='__main__':
    Entrez.email = 'ekolomenski@gmail.com';
    Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';
    mm9 = EnsemblRelease(67, species='mouse');
    hg19 = EnsemblRelease(102, species='human');
    GRCm38 = EnsemblRelease(102, species='mouse');
    dict_genomas = {'mm9':mm9, 'mouse':mm9, 'hg19':hg19, 'human':hg19, 'GRCm38':GRCm38, 'mouse102':GRCm38};
    output_dump.append('');