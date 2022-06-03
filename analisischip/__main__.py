import logging
from analisischip import seq_data

'''
Archivo donde ejecuto las pruebas
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
    output_dump.append('');