# Generales
import os
import time
import logging
# Analisis de secuencias y genomas
from Bio import Entrez, SeqIO


###################################### CLASES #####################################


class seq_data(object):
    '''
Clase para cargar datos de secuencia en un genoma dado
Almacena secuencias en formato chr_n, pos_ini, pos_end
Descarga archivos .fasta con los cromosomas necesarios para las secuencias usadas
    '''
    def __init__(self, genome_name, genome_element, path_fasta=''):
        # M_seq almacena todos los rangos de secuencias
        # Rangos almacenados en formato chr_n, pos_ini, pos_end
        # Secuencia almacenada en archivos en path a definir

        # dict_range registra los rangos por chr_n
        self.dict_range = {};
        # Registro el genoma y su nombre
        self.genome_name = genome_name;
        self.genome = genome_element;
        self.path_fasta = path_fasta;
        if path_fasta == '':
            logging.warning('No se definio path_fasta. Se buscan y descargan archivos de secuencia en directorio actual.');
        # Cargo dict_chrid para pasar de chr_n a chr_id
        self._cargar_dict_chrid(genome_name);
        return None


    ## FUNCIONES:
    # X cargar_rango(chr_n, pos_ini, pos_end): Carga un rango de secuencias en self.dict_range (y hace muchos chequeos)
    # X _chr_check(chr_n): Al cargar, ver si el cromosoma esta presente en self.dict_range o en carpeta path_fasta
    # X _chr_file_check(chr_n): Busca el archivo en carpeta path_fasta y usa _download_chr(chr_n) si no esta
    # X _download_chr(chr_n, retries): Ver consulta_secuencia_chr() en 14-PruebaDescargarChr.py
    # X _buscar_chrid(chr_n): Consigue chrID en base a chr_n, usando self.dict_chrid cargado en __init__
    # X _cargar_dict_chrid(genome_name): Carga self.dict_chrid para pasar de chr_n a chr_id
    # X _agregar_chrid(chr_id, chr_n): Revisa que chr_n no este en self.dict_chrid y agrega chr_id si no esta
    # X _consulta_entrez_chr(chr_id): Consigue el elemento correspondiente al cromosoma con SeqIO
    # ? _check_overlap(chr_n, pos_ini, pos_end): Revisa que no haya overlap (ver funciones armadas en este archivo)
    # - cargar_bed(archivo): Carga todos los rangos en un archivo de output de ChIP-seq
    # - cargar_promotores_genoma(rango): Carga todos los rangos alrededor de promotores de genes
    # - Con funcion cargar_rango() funcional, hacer cargar_bed() y cargar_promotores_genoma(rango)


    def _agregar_chrid(self, chr_id, chr_n):
        # Revisa que chr_n este en self.dict_chrid.keys() y agrega chr_id a self.dict_chrid[chr_n]
        # Si ya se encuentra chr_n en self.dict_chrid.keys() y chr_id es distinto, tira warning y actualiza chr_id

        # Si chr_n esta en self.dict_chrid.keys(), revisa que tenga el mismo chr_id
        if chr_n in self.dict_chrid.keys():
            if self.dict_chrid[chr_n] == chr_id:
                print('Id "' + chr_id + '" ya asociado con "' + chr_n + '".')
            else:
                logging.warning('"' + chr_n + '" asociado con id "' + self.dict_chrid[chr_n] + '". Se reemplaza por id "' + chr_id + '".')
        # Si chr_n no esta en self.dict_chrid.keys(), se agrega directamente
        else:
            self.dict_chrid[chr_n] = chr_id;
            print('Dict chrid actualizado correctamente con id "' + chr_id + '" para "' + chr_n + '".')
        return self


    def _buscar_chrid(self, chr_n):
        # Busca el id de nucleotide para chr_n en self.genome

        # Veo que chr_n este en self.dict_chrid.keys(), por las dudas
        if chr_n in self.dict_chrid.keys():
            ret = self.dict_chrid[chr_n];
        # Tira error si no esta y devuelve string vacio
        else:
            logging.error('No se pudo encontrar ' + chr_n + ' en dict_chrid. Se puede agregar manualmente con self._agregar_chrid(chr_id, chr_n).')
            ret = '';
        return ret


    def _cargar_dict_chrid(self, genome_name):
        # Inicializa self.dict_chrid para pasar de chr_n a chr_id
        # Redefine self.genome_name para estandarizar outputs
        self.dict_chrid = {};
        if genome_name.lower() == 'hg19' or genome_name.lower() == 'human':
            self.dict_chrid = {'chr1':'NC_000001.11', 'chr2':'NC_000002.12', 'chr3':'NC_000003.12', 'chr4':'NC_000004.12',
                               'chr5':'NC_000005.10', 'chr6':'NC_000006.12', 'chr7':'NC_000007.14', 'chr8':'NC_000008.11',
                               'chr9':'NC_000009.12', 'chr10':'NC_000010.11', 'chr11':'NC_000011.10', 'chr12':'NC_000012.12',
                               'chr13':'NC_000013.11', 'chr14':'NC_000014.9', 'chr15':'NC_000015.10', 'chr16':'NC_000016.10',
                               'chr17':'NC_000017.11', 'chr18':'NC_000018.10', 'chr19':'NC_000019.10', 'chr20':'NC_000020.11',
                               'chr21':'NC_000021.9', 'chr22':'NC_000022.11', 'chrX':'NC_000023.11', 'chrY':'NC_000024.10',
                               'chrM':'NC_012920.1', 'chrMT':'NC_012920.1'};
            self.genome_name = 'hg19';
        elif genome_name.lower() == 'mm9' or genome_name.lower() == 'mouse':
            self.dict_chrid = {'chr1':'NC_000067.5', 'chr2':'NC_000068.6', 'chr3':'NC_000069.5', 'chr4':'NC_000070.5',
                               'chr5':'NC_000071.5', 'chr6':'NC_000072.5', 'chr7':'NC_000073.5', 'chr8':'NC_000074.5',
                               'chr9':'NC_000075.5', 'chr10':'NC_000076.5', 'chr11':'NC_000077.5', 'chr12':'NC_000078.5',
                               'chr13':'NC_000079.5', 'chr14':'NC_000080.5', 'chr15':'NC_000081.5', 'chr16':'NC_000082.5',
                               'chr17':'NC_000083.5', 'chr18':'NC_000084.5', 'chr19':'NC_000085.5', 'chrX':'NC_000086.6',
                               'chrY':'NC_000087.6', 'chrM':'NC_005089.1', 'chrMT':'NC_005089.1'};
            self.genome_name = 'mm9';
        elif genome_name.lower() == 'hg38':
            self.dict_chrid = {'chr1':'CM000994.3', 'chr10':'CM001003.3', 'chr11':'CM001004.3', 'chr12':'CM001005.3',
                               'chr13':'CM001006.3', 'chr14':'CM001007.3', 'chr15':'CM001008.3', 'chr16':'CM001009.3',
                               'chr17':'CM001010.3', 'chr18':'CM001011.3', 'chr19':'CM001012.3', 'chr2':'CM000995.3',
                               'chr3':'CM000996.3', 'chr4':'CM000997.3', 'chr5':'CM000998.3', 'chr6':'CM000999.3',
                               'chr7':'CM001000.3', 'chr8':'CM001001.3', 'chr9':'CM001002.3', 'chrMT':'AY172335.1',
                               'chrX':'CM001013.3', 'chrY':'CM001014.3', 'chrM':'AY172335.1'};
            self.genome_name = 'hg38';
        else:
            logging.warning('No se encontro el genoma ' + genome_name + ', se asume genoma hg19 (Homo sapiens).');
            self.dict_chrid = {'chr1':'NC_000001.11', 'chr2':'NC_000002.12', 'chr3':'NC_000003.12', 'chr4':'NC_000004.12',
                               'chr5':'NC_000005.10', 'chr6':'NC_000006.12', 'chr7':'NC_000007.14', 'chr8':'NC_000008.11',
                               'chr9':'NC_000009.12', 'chr10':'NC_000010.11', 'chr11':'NC_000011.10', 'chr12':'NC_000012.12',
                               'chr13':'NC_000013.11', 'chr14':'NC_000014.9', 'chr15':'NC_000015.10', 'chr16':'NC_000016.10',
                               'chr17':'NC_000017.11', 'chr18':'NC_000018.10', 'chr19':'NC_000019.10', 'chr20':'NC_000020.11',
                               'chr21':'NC_000021.9', 'chr22':'NC_000022.11', 'chrX':'NC_000023.11', 'chrY':'NC_000024.10',
                               'chrM':'NC_012920.1', 'chrMT':'NC_012920.1'};
            self.genome_name = 'hg19';
        return self


    def _chr_check(self, chr_n):
        # Revisa si chr_n esta presente en self.dict_range o en carpeta self.path_fasta
        # Si chr_n esta presente o se puede agregar, devuelve True
        # Si algo falla, devuelve False

        # Inicializo el booleano que se devuelve
        ret = False;
        # Si chr_n esta en self.dict_range.keys(), no hago nada y devuelvo True
        if chr_n in self.dict_range.keys():
            ret = True;
        # Si chr_n no esta en self.dict_range.keys(), reviso carpeta self.path_fasta con self._chr_file_check()
        elif self._chr_file_check(chr_n):
            # Una vez me aseguro que el archivo esta presente, agrego chr_n a self.dict_range
            self.dict_range[chr_n] = [];
            ret = True;
        return ret


    def _chr_file_check(self, chr_n):
        # Busca el archivo correspondiente a chr_n en carpeta path_fasta 
        # Devuelve True si el archivo esta presente o si se descarga correctamente
        # Devuelve False si el archivo no esta presente y no se pudo descargar por una u otra razon

        # Defino el nombre y direccion del archivo a buscar
        nom_arch = self.genome_name + '_' + chr_n + '.fasta';
        if self.path_fasta == '':
            dir_arch = '.\\' + nom_arch;
        else:
            dir_arch = os.path.join(self.path_fasta, nom_arch);
        # Busco el archivo self.genome_name + '_' + chr_n + '.fasta' en self.path_fasta
        try:
            # Si esta presente, no hago nada
            F = open(dir_arch, 'r');
            F.close();
            ret = True;
        except:
            logging.warning('Archivo ' + nom_arch + ' no encontrado en la carpeta dada.');
            # Si chr_n no esta presente, defino chr_id en base a chr_n para descargarlo
            chr_id = self._buscar_chrid(chr_n);
            # Si no encuentro chr_id, tiro error y devuelvo False
            if chr_id == '':
                logging.error('Cromosoma ' + chr_n + ' no encontrado en lista de IDs para buscar en nucleotide.');
                ret = False;
            # Si encuentro chr_id, trato de bajarlo con self._download_chr()
            else:
                ret = self._download_chr(chr_n);
        return ret


    def _consulta_entrez_chr(self, chr_id):
        # Devuelve objeto record de Entrez conteniendo la secuencia completa de un cromosoma en base al ID de base de datos nucleotide

        # Busco el ID en la base de datos nucleotide
        handle_nuc_id = Entrez.esearch(db='nucleotide', term=chr_id, retmax=1000);
        record_nuc_id = Entrez.read(handle_nuc_id);
        nuc_id = record_nuc_id['IdList'];
        # Reviso que se haya encontrado un solo ID, sino tiro error y trabajo solo con el primero
        if len(nuc_id) > 1:
            logging.warning('Mas de un resultado para el id "' + chr_id + '". Se usa solo el primero.');
            print('Lista completa: \n' + str(nuc_id))
        # Uso el nuc_id para conseguir el record del cromosoma entero en nucleotide
        handle = Entrez.efetch(db='nucleotide', rettype='gbwithparts', retmode='text', id=nuc_id[0]);
        record = SeqIO.read(handle, 'genbank');
        return record


    def _download_chr(self, chr_n, retries=10):
        # Descarga el archivo .fasta correspondiente a chr_n para self.genome_name
        # Devuelve True si puede descargar el archivo
        # Devuelve False si falla la descarga

        # Inicializo el booleano que se devuelve
        ret = False;
        # Chequeo para retries < 1
        if retries < 1:
            retries = 1;
        # Consigo chr_id en base a chr_n
        chr_id = self._buscar_chrid(chr_n);
        # Inicializo las variables del ciclo while
        encontrado = False;
        tries = 0;
        # Repito maximo retries veces o hasta que se encuentre
        while not encontrado and tries < retries:
            try:
                # Uso _consulta_entrez_chr(chr_id) para conseguir el objeto record con los datos para el .fasta
                record_chr = self._consulta_entrez_chr(chr_id);
                # Si record_chr tiene largo mayor a 0, lo guardo en archivo .fasta dentro de self.path_fasta
                if len(record_chr) > 0:
                    print('Objeto record de Entrez obtenido para id "' + chr_id + '". Guardando en archivo .fasta.')
                    # Defino el nombre del archivo
                    nom_seqio = self.genome_name + '_' + chr_n + '.fasta';
                    # Defino el path completo para el archivo
                    if self.path_fasta == '':
                        path_seqio = '.\\' + nom_seqio;
                    else:
                        path_seqio = os.path.join(self.path_fasta, nom_seqio);
                    # Uso SeqIO para guardarlo
                    SeqIO.write(record_chr, path_seqio, 'fasta');
                    print('Archivo .fasta creado correctamente.')
                    # Una vez guardado, registro encontrado = True y ret = True
                    encontrado = True;
                    ret = True;
                # Si record_chr no encuentra nada (largo=0), tiro warning y devuelvo false
                else:
                    logging.error('No se pudo encontrar secuencia para "' + chr_id + '".')
                    # No registro encontrado = True por si hay un error en descarga y un retry lo puede solucionar
            except:
                tries = tries + 1;
                print('Fallo intento ' + str(tries) + '.')
        return ret


    def cargar_rango(self, chr_n, pos_ini, pos_end):
        # Carga un rango de secuencias en self.dict_range (y hace muchos chequeos)
        # No revisa si hay superposicion

        # Uso self._chr_check(chr_n) para revisar que el cromosoma este presente
        chr_presente = self._chr_check(str(chr_n));
        # Si esta presente, trato de agregar pos_ini y pos_end
        if chr_presente:
            # Es posible que sea necesario revisar superposicion, pero inicialmente no lo hago
            # Creo el rango revisando que sean int y pos_ini sea menor que pos_end
            rango_cargado = (min(int(pos_ini), int(pos_end)), max(int(pos_ini), int(pos_end)));
            self.dict_range[chr_n].append(rango_cargado);
        # Si no pude agregar chr_n, tiro error
        else:
            logging.error('Rango no agregado por no poder parsear "' + str(chr_n) + '"')
        return self


#################################### FUNCIONES ####################################


def _main_test():
    # Funcion para probar funciones en ejecucion del archivo

    L_out = [];
    
    return L_out


##################################### PRUEBAS #####################################


### Funciones para pruebas

def ConsultaSecuencia(id_chr, seq_start, seq_finish, strand=1, sleep_time=60):
    # Devuelve una secuencia dado un ID de cromosoma (incluye info de especie) y posicion inicial/final
    time.sleep(0.1);
    rec_seq = '';
    try:
        handle = Entrez.efetch(db='nucleotide', id=id_chr, rettype='fasta',
                               strand=strand, seq_start=seq_start, seq_stop=seq_finish);
        record = SeqIO.read(handle, 'fasta');
        handle.close();
        rec_seq = record.seq;
    except:
        logging.warning('Exception raised for chr ' + str(id_chr) + ' between positions ' + 
                        str(seq_start) + ' and ' + str(seq_finish) + '.');
        time.sleep(sleep_time);
        try:
            handle = Entrez.efetch(db='nucleotide', id=id_chr, rettype='fasta',
                                   strand=strand, seq_start=seq_start, seq_stop=seq_finish);
            record = SeqIO.read(handle, 'fasta');
            handle.close();
            rec_seq = record.seq;
        except:
            logging.error('Retry failed. Returning empty string.');
    return rec_seq


def IDchr(chromosome,genome='hg19'):
    # Devuelve el ID de cromosoma para ConsultaSecuencia dado el numero de cromosoma y el genoma correspondiente
    # Programado a mano, funciona solo con hg19 (humano) y mm9 (raton)
    ret = '';
    b = True;
    if genome.lower() == 'hg19' or genome.lower() == 'human': # GRCh38.p13
        dict_IDchr = {'1':'NC_000001.11', '2':'NC_000002.12', '3':'NC_000003.12', '4':'NC_000004.12',
                      '5':'NC_000005.10', '6':'NC_000006.12', '7':'NC_000007.14', '8':'NC_000008.11',
                      '9':'NC_000009.12', '10':'NC_000010.11', '11':'NC_000011.10', '12':'NC_000012.12',
                      '13':'NC_000013.11', '14':'NC_000014.9', '15':'NC_000015.10', '16':'NC_000016.10',
                      '17':'NC_000017.11', '18':'NC_000018.10', '19':'NC_000019.10', '20':'NC_000020.11',
                      '21':'NC_000021.9', '22':'NC_000022.11', 'X':'NC_000023.11', 'Y':'NC_000024.10',
                      'M':'NC_012920.1', 'MT':'NC_012920.1'};
    elif genome.lower() == 'mm9' or genome.lower() == 'mouse': # MGSCv37
        dict_IDchr = {'1':'NC_000067.5', '2':'NC_000068.6', '3':'NC_000069.5', '4':'NC_000070.5',
                      '5':'NC_000071.5', '6':'NC_000072.5', '7':'NC_000073.5', '8':'NC_000074.5',
                      '9':'NC_000075.5', '10':'NC_000076.5', '11':'NC_000077.5', '12':'NC_000078.5',
                      '13':'NC_000079.5', '14':'NC_000080.5', '15':'NC_000081.5', '16':'NC_000082.5',
                      '17':'NC_000083.5', '18':'NC_000084.5', '19':'NC_000085.5', 'X':'NC_000086.6',
                      'Y':'NC_000087.6', 'M':'NC_005089.1', 'MT':'NC_005089.1'};
    elif genome.lower() == 'mouse102' or genome.lower() == 'grcm38':
        dict_IDchr = {'1':'CM000994.3', '10':'CM001003.3', '11':'CM001004.3', '12':'CM001005.3',
                      '13':'CM001006.3', '14':'CM001007.3', '15':'CM001008.3', '16':'CM001009.3',
                      '17':'CM001010.3', '18':'CM001011.3', '19':'CM001012.3', '2':'CM000995.3',
                      '3':'CM000996.3', '4':'CM000997.3', '5':'CM000998.3', '6':'CM000999.3',
                      '7':'CM001000.3', '8':'CM001001.3', '9':'CM001002.3', 'MT':'AY172335.1',
                      'X':'CM001013.3', 'Y':'CM001014.3', 'M':'AY172335.1'};
    else:
        logging.error('No se pudo encontrar genoma ' + str(genome));
        b = False;
        dict_IDchr = {};

    if str(chromosome).upper() in dict_IDchr.keys():
        ret = dict_IDchr[str(chromosome).upper()];
    elif b:
        logging.error('No se pudo encontrar cromosoma ' + str(chromosome));
    return ret


###################################################################################
######################################## OLD ######################################
###################################################################################


###################################### CLASES #####################################
'''
    def _revisar_overlap_dict(self, chr_n, pos_ini, pos_end, seq):
        # Revisa que no haya overlap en self.dict_range[chr_n] para pos_ini, pos_end
        # Primero carga el rango [pos_ini, pos_end] en orden
        # Despues revisa que no haya overlap con secuencia anterior y/o siguiente
        # Si hay overlap, extiende la secuencia correspondiente, une dos rangos o ignora el rango dado

        ### Display
        print('Revisando overlap en ' + str(chr_n) + ' para rango ' + str(pos_ini) + ' a ' + str(pos_end))

        # Booleano para saber si se encontro la posicion para el rango
        posicion_encontrada = False;
        # Contador para recorrer self.dict_range[chr_n]
        i = 0;
        largo_dict_chr_n = len(self.dict_range[chr_n]);
        # Booleano para chequeos de overlap
        ultima_posicion = False;
        # Ciclo para recorrer self.dict_range[chr_n] y agregar [pos_ini, pos_end] en el lugar correcto
        while (not posicion_encontrada) and i < largo_dict_chr_n:
            curr_range = self.dict_range[chr_n][i];
            # Si pos_ini es menor que el rango en posicion i, se aumenta i en 1
            if pos_ini > curr_range[0]:
                i = i+1;
                # Si se llega a la ultima posicion, se agrega [pos_ini, pos_end] al final
                if i == largo_dict_chr_n:
                    self.dict_range[chr_n].append([pos_ini, pos_end]);
                    ultima_posicion = True;
            # En el momento que se encuentre un rango con pos_ini mayor al dado, se carga el rango
            else:
                posicion_encontrada = True;
                # Se agrega [pos_ini, pos_end] en la posicion i, moviendo el rango que estaba en esa posicion para arriba
                self.dict_range[chr_n] = self.dict_range[chr_n][:i] + [[pos_ini, pos_end]] + self.dict_range[chr_n][i:];

        # Reviso posiciones adyacentes a self.dict_range[chr_n][i] (+1 y -1)
        overlap_antes = False;
        overlap_despues = False;
        nueva_pos_ini = None;
        nueva_pos_end = None;
        # Si el rango agregado esta en la ultima posicion, no hay rango siguiente y no se hace este chequeo
        if not ultima_posicion:
            ##print('Revisando overlap_despues')
            ##print('self.dict_range[chr_n][i+1][0] para i=' + str(i)+ ' es ' + str(self.dict_range[chr_n][i+1][0]))
            ##print('Rango: (' +str(pos_ini) + ',' + str(pos_end) + ')')
            # Si self.dict_range[chr_n][i+1][0] es menor a pos_end, hay overlap con la secuencia siguiente
            if self.dict_range[chr_n][i+1][0] < pos_end:
                overlap_despues = True;
                nueva_pos_ini = pos_ini;
                nueva_pos_end = max(self.dict_range[chr_n][i+1][1], pos_end);
        # Si i es 0, no hay rango anterior y no se hace este chequeo
        if i > 0:
            ##print('Revisando overlap_antes')
            ##print('self.dict_range[chr_n][i-1][1] para i=' + str(i)+ ' es ' + str(self.dict_range[chr_n][i-1][1]) + ')')
            ##print('Rango: (' +str(pos_ini) + ',' + str(pos_end))
            # Si self.dict_range[chr_n][i-1][1] es mayor a pos_ini, hay overlap con la secuencia anterior
            if self.dict_range[chr_n][i-1][1] > pos_ini:
                overlap_antes = True;
                nueva_pos_ini = self.dict_range[chr_n][i-1][0];
                nueva_pos_end = max(self.dict_range[chr_n][i-1][1], pos_end);
        
        # Si hay overlap antes y despues, uno los tres rangos en uno
        if overlap_antes and overlap_despues:
            nueva_pos_ini = self.dict_range[chr_n][i-1][0];
            nueva_pos_end = max(self.dict_range[chr_n][i-1][1], pos_end, self.dict_range[chr_n][i+1][1]);
            self.dict_range[chr_n] = self.dict_range[chr_n][:i-1] + [[nueva_pos_ini, nueva_pos_end]] + self.dict_range[chr_n][i+2:];
        # Si solo hay overlap antes o despues, uno dos rangos
        elif overlap_antes:
            self.dict_range[chr_n] = self.dict_range[chr_n][:i-1] + [[nueva_pos_ini, nueva_pos_end]] + self.dict_range[chr_n][i+1:];
        elif overlap_despues:
            self.dict_range[chr_n] = self.dict_range[chr_n][:i] + [[nueva_pos_ini, nueva_pos_end]] + self.dict_range[chr_n][i+2:];
        # Si no hay overlap ni antes ni despues, el rango queda agregado

        ### Display
        ##print('Overlap antes: ' + str(overlap_antes))
        ##print('Overlap despues: ' + str(overlap_despues))
        return self


    def actualizar_M_seq(self, buscar_seq):
        # Actualiza M_seq en base a los rangos registrados en dict_range

        ########### FALTA
        return self


    def agregar_secuencia(self, chr_n, pos_ini, pos_end, seq=''):
        # Agrega una secuencia a self.dict_range
        # Revisa que no haya overlap con ninguna secuencia en self.dict_range
        # Si hay overlap, hace chequeo y merge
        # Si no hay overlap, agrega el rango al diccionario

        ### Display
        print('Empezando parseo de rango ' + str(chr_n) + ' : ' + str(pos_ini) + ' - ' + str(pos_end))

        # Parseo el rango para que pos_ini sea el menor numero y pos_end, el mayor
        seq_range = [min(int(pos_ini), int(pos_end)), max(int(pos_ini), int(pos_end))];
        # Primero reviso si chr_n esta en self.dict_range.keys()
        if chr_n in self.dict_range.keys():
            # Si chr_n ya esta en self.dict_range.keys(), tengo que chequear que no haya overlap
            self._revisar_overlap_dict(chr_n, seq_range[0], seq_range[1], seq);
        else:
            print('Rango agregado a dict_range directamente')
            # Si chr_n no esta en self.dict_range.keys(), agrego la key y el rango
            self.dict_range[chr_n] = [];
            self.dict_range[chr_n].append(seq_range);

        return self


    def _buscar_secuencia(self, M_seq_id):
        # Busca la secuencia para self.M_seq[M_seq_id]
        # Si seq == '', la carga directamente
        # Si seq == secuencia buscada, da OK
        # Si seq != secuencia buscada, tira error y cambia la secuencia en self.M_seq
        # Devuelve un string que se printea con _print_progress() si el ciclo quiere

        # Inicializo la variable que se devuelve
        p_ret = '';
        # Defino el sitio y la secuencia
        curr_site = self.M_seq[M_seq_id];
        seq_site = curr_site[-1];
        # Busco el id del cromosoma
        chr_id = IDchr(curr_site[0][3:],genome=self.genome_name);
        # Busco la secuencia dada por pos_ini y pos_end
        curr_seq = str(ConsultaSecuencia(chr_id, curr_site[1], curr_site[2]));

        if seq_site == '':
            self.M_seq[M_seq_id][-1] = str(curr_seq);
            p_ret = 'Sitio vacio encontrado y cargado';
        elif seq_site != curr_seq:
            p_ret = 'ERROR: Sitio encontrado distinto al cargado';
        else:
            p_ret = 'OK';
        return p_ret


    def _cargar_rango(self, chr_n, pos_ini, pos_end, seq=''):
        # Agrega chr_n, pos_ini, pos_end y seq directamente a self.M_seq
        # Funcion a fuerza bruta para no usar por si sola
        self.M_seq.append([chr_n, pos_ini, pos_end, seq]);
        return self


    def _print_progress(self, texto):
        # Funcion para hacer prints que se sobreescriban entre si
        # Cambio end para usar en consola o en Python IDLE
        # end='\r' para consola; end='\n' para Python IDLE
        print(texto, end='\r');
        return self


    def agregar_rango(self, chr_n, pos_ini, pos_end, seq=''):
        # Agrega rango en formato chr_n, pos_ini, pos_end a self.M_seq
        # Funcion que revisa si el rango ya existe en self.M_seq o si hay overlap
        # Si el rango ya esta representado, no hace nada
        # Si el rango no esta representado ni toca con otros, lo agrega directamente
        # Si el rango tiene overlap con otro rango, los junta y transforma en un rango mayor

        # Primero veo si hay overlap
        overlap_check = self.buscar_overlap(chr_n, pos_ini, pos_end); ####### FALTA HACER ESTO #######

        ##################### HACER #####################
        return self


    def append_M_seq(self, loaded_M):
        # Funcion que revisa una matriz cargada y la agrega a self.M_seq

        # Reviso cada rango en loaded_M
        for i in range(len(loaded_M)):
            curr_range = loaded_M[i];
            if len(curr_range)>3:
                # Dejo que append_range() decida si se agrega, se solapa o se ignora
                self.agregar_rango(curr_range[0], curr_range[1], curr_range[2], seq=curr_range[3]); ####### FALTA HACER ESTO #######
            else:
                # Dejo que append_range() decida si se agrega, se solapa o se ignora
                self.agregar_rango(curr_range[0], curr_range[1], curr_range[2], seq=''); ####### FALTA HACER ESTO #######
        return self


    def buscar_overlap(self, chr_n, pos_ini, pos_end):
        # Busca si self.M_seq tiene segmentos que se solapen con chr_n, pos_ini, pos_end
        # Devuelve 0 si el rango no esta en self.M_seq
        # Devuelve 1 si hay overlap o si el rango toca a otro rango en self.M_seq
        # Devuelve 2 si el rango se encuentra dentro de self.M_seq
        
        # Inicializo los booleanos que definen que pasa
        # M_seq_overlap registra todos los rangos que se superpongan con el buscado
        M_seq_overlap = [];
        # Recorro self.M_seq
        for i in range(len(self.M_seq)):
            # Primero selecciono los rangos de self.M_seq en chr_n
            if self.M_seq[i][0] == chr_n:
                ########### HACER CHEQUEO DE OVERLAP ACA ###########
                M_seq_overlap.append(self.M_seq[i][:]);
        
        ##################### HACER #####################
        ## FALTA:
        ## Ver cuantos elementos hay en M_seq_overlap
        ## Si hay uno, ver si es overlap completo o parcial
        ## Si hay cero, agregar el rango dado
        ## Si hay mas de uno, ver si tienen overlap entre ellos primero
        ## Si no tienen overlap, hacer chequeo fuerza bruta
        return self


    def cargar_archivo(self, nombre_in, append_mode=False, ext='.csv', sep=';'):
        # Funcion para cargar todos los rangos y secuencias desde un archivo
        # Optimizado para funcionar con el output de self.guardar_archivo()
        # append_mode define si se agregan los sitios del archivo o si se borra todo antes de cargar

        # Extraigo las secuencias de nombre_in en loaded_M_seq
        loaded_M_seq = abrir_archivo(nombre_in, ext, sep);

        # Si append_mode es True, se compara self.M_seq con loaded_M_seq
        if append_mode:
            self.append_M_seq(loaded_M_seq); ############# FALTA HACER #############
        # Si append_mode es False, se reinicia self.M_seq con loaded_M_seq
        else:
            self.M_seq = loaded_M_seq;
        return self


    def leer_bed(self, nom_arch, append_mode=True, path_arch='', ext='.bed', sep='\t'):
        # Funcion para leer archivos .bed y guardar los rangos de peaks en self.M_seq
        # Asume que chr_n, pos_ini y pos_end estan en las primeras 3 filas
        ### Puedo agregar L_data_id para definir las columnas con chr_n, pos_ini y pos_end

        # Defino dir_arch con nom_arch y path_arch
        if path_arch != '':
            dir_arch = os.path.join(path_arch, nom_arch);
        else:
            dir_arch = nom_arch;
        # Extraigo las secuencias de dir_arch en loaded_M_seq
        loaded_M_seq = abrir_arch_dir(dir_arch, nom_arch, ext, sep);

        # Selecciono las primeras 3 columnas de cada fila
        for i in range(len(loaded_M_seq)):
            loaded_M_seq[i] = loaded_M_seq[i][:3];

        # Si append_mode es True, se compara self.M_seq con loaded_M_seq
        if append_mode:
            self.append_M_seq(loaded_M_seq); ############# FALTA HACER #############
        # Si append_mode es False, se reinicia self.M_seq con loaded_M_seq
        else:
            self.M_seq = loaded_M_seq;
        return self


    def guardar_archivo(self, nombre_out, ext='.csv', sep=';'):
        # Funcion para guardar todos los rangos y secuencias en un .csv
        
        # Creo y abro el archivo de output
        with open(nombre_out + ext, 'w') as F_out:
            print('Archivo ' + str(nombre_out) + str(ext) + ' creado.');
        with open(nombre_out + ext, 'a') as F_out:
            # Recorro self.M_seq
            for i in range(len(self.M_seq)):
                curr_SU = self.M_seq[i];
                # Inicializo el texto que se guarda en el archivo
                curr_r = '';
                # Cargo cada elemento de curr_SU separado por sep
                for j in curr_SU:
                    curr_r = curr_r + str(j) + str(sep);
                # Elimino el ultimo sep y agrego el fin de linea
                curr_r = curr_r.rstrip(sep) + '\n';
                # Guardo el texto en F_out
                F_out.write(str(curr_r));
        return self


    def revisar_mult(self, carga_vacias=False, L_ids=[]):
        # Revisa si self.M_seq[id] contiene seq correspondiente al rango
        # Consulta por secuencias a Ensembl (lento para muchas consultas)
        # Si L_ids es vacio, revisa TODO (puede llegar a ser MUY lento)
        # Puede usarse para cargar secuencias vacias

        # Primero reviso si L_ids esta vacio
        if len(L_ids) == 0:
            cont_total = len(self.M_seq);
            # Si esta vacio, armo un L_ids con todos los ids
            L_ids_usado = list(range(cont_total));
        else:
            cont_total = len(L_ids);
            # Sino uso L_ids directamente
            L_ids_usado = L_ids;

        # Variables para display
        seq_ok = 0;
        seq_error = 0;
        seq_cargada = 0;
        seq_vacias = 0;
        # Recorro L_ids_usado
        for i in L_ids_usado:
            # Texto para display
            p_text = '';
            # Defino curr_sitio y curr_seq
            curr_site = self.M_seq[i];
            curr_seq = curr_site[-1];
            # Si la secuencia esta vacia, defino que se hace con carga_vacias e ignora_vacias
            if curr_seq == '':
                # Si carga_vacias == True, uso self._buscar_secuencia() 
                if carga_vacias:
                    p_ret = self._buscar_secuencia(i);
                    seq_cargada = seq_cargada + 1;
                # Sino solo registro seq_vacias
                else:
                    seq_vacias = seq_vacias + 1;
            # Si la secuencia no esta vacia y es del largo esperado por pos_ini y pos_end
            # Paso directamente a self._buscar_secuencia()
            elif len(curr_seq) == curr_site[2]-curr_site[1]+1:
                p_ret = self._buscar_secuencia(i);
                if p_ret == 'OK':
                    seq_ok = seq_ok + 1;
                else:
                    seq_error = seq_error + 1;
            # Si la secuencia es de distinto largo del esperado
            # Paso a self._buscar_secuencia() y agrego directamente a seq_error
            else:
                p_ret = self._buscar_secuencia(i);
                p_ret = p_ret + ' y de distinto largo';
                seq_error = seq_error + 1;
            ### Display
            cont = seq_ok+seq_error+seq_cargada+seq_vacias;
            p_text = 'Revisadas ' + str(cont) + ' de ' + str(len(L_ids_usado)) + ' secuencias.';
            self._print_progress(p_text);
            # Reviso si hubo errores que mostrar
            if p_ret[5:] == 'ERROR':
                logging.error(p_ret);
        p_final = 'Carga terminada. Secuencias OK: ' + str(seq_ok) + '. Secuencias con errores: ' + str(seq_error) + '.'
        if carga_vacias:
            p_final = p_final + ' Secuencias cargadas: ' + str(seq_cargada) + '.';
        else:
            p_final = p_final + ' Secuencias vacias: ' + str(seq_vacias) + '.';
        print(p_final);
        return self
'''

#################################### FUNCIONES ####################################


def abrir_arch_dir(dir_arch, nom_arch, ext, sep_arch):
    # Abre un archivo en una carpeta diferente a la actual, similar a abrir_archivo()
    # Devuelve las filas divididas en columnas como matriz

    # Creo la matriz que se devuelve
    M_out = [];
    # Abro el archivo
    with open(dir_arch + ext, 'r') as F:
        print('Archivo ' + str(nom_arch) + str(ext) + ' abierto.');
        # Reviso cada linea de F
        for curr_line in F:
            # Paso cada linea a formato de lista
            L = curr_line.rstrip().split(sep_arch);
            # Guardo L en M_out
            M_out.append(L[:]);
    return M_out


def abrir_archivo(nom_arch, ext, sep_arch):
    # Abre un archivo y devuelve las filas divididas en columnas como matriz

    # Creo la matriz que se devuelve
    M_out = [];
    # Abro el archivo
    with open(nom_arch + ext, 'r') as F:
        print('Archivo ' + str(nom_arch) + str(ext) + ' abierto.');
        # Reviso cada linea de F
        for curr_line in F:
            # Paso cada linea a formato de lista
            L = curr_line.rstrip().split(sep_arch);
            # Guardo L en M_out
            M_out.append(L[:]);
    return M_out


def complemento(N,adn=True):
    # Devuelve el complemento de un nucleotido en adn o arn
    dict_adn = {'T':'A','U':'A','A':'T','C':'G','G':'C','N':'N'};
    dict_arn = {'T':'A','U':'A','A':'U','C':'G','G':'C','N':'N'};

    if not (N in dict_adn.keys()):
        logging.warning('Nucleotido "' + str(N) + '" no interpretado. Se devuelve N.');
        ret = 'N';
    elif adn:
        ret = dict_adn[N];
    else:
        ret = dict_arn[N];
    return(ret)


def complemento_secuencia(seq, adn=True):
    # Devuelve el complemento de una secuencia de adn o arn
    # La secuencia del complemento se devuelve en la misma orientacion (al reves que la referencia)
    ret_seq = '';
    for i in range(len(seq)):
        ret_seq = ret_seq + complemento(seq[-i-1],adn=adn);
    return(ret_seq)


def max_range(M_num):
    # Recibe una matriz que sea una lista de listas de numeros
    # Devuelve el rango que cubra todos los numeros en M_num

    # Inicializo los valores inicial y final
    num_min = '';
    num_max = '';
    # Recorro M_num y busco el mayor y el menor numero
    for i in M_num:
        if num_min == '':
            num_min = min(i);
        else:
            num_min = min(min(i), num_min);
        if num_max == '':
            num_max = max(i);
        else:
            num_max = max(max(i),num_max);
    return [num_min, num_max]


def merge_overlapping_intervals(L_intervals):
    # Recibe una lista de intervalos y devuelve una lista de intervalos sin overlap
    # Une los intervalos que se superpongan

    # Ordeno la lista de rangos por su primer elemento
    L_in = sorted(L_intervals, key=lambda x: x[0]);
    # Inicializo la variable que se devuelve
    L_out = [];

    # Solo agrego algo a L_out si L_in tiene mas de un elemento
    if len(L_in)>0:
        L_out.append(L_in[0][:]);
        # Recorro L_in a partir del segundo elemento
        for i in range(1,len(L_in)):
            # Defino el proximo rango 
            curr_range = L_in[i];
            # Agarro el ultimo elemento de L_out (para hacer merge de ser necesario)
            pop_range = L_out.pop();
            # Si hay overlap, hago merge
            if range_overlap(pop_range, curr_range):
                # merged_range siempre empieza con pop_range[0] porque esta ordenado de esa manera
                merged_range = [pop_range[0], max(pop_range[1], curr_range[1])];
                L_out.append(merged_range[:]);
            # Si no hay overlap, devuelvo pop_range a L_out y agrego curr_range al final
            else:
                L_out.append(pop_range[:]);
                L_out.append(curr_range[:]);
    return L_out


def range_overlap(range1, range2):
    # Revisa si dos rangos tienen overlap
    # range1[0] es el numero mas bajo dado
    return range1[1] >= range2[0]


###################################################################################
####################################### MAIN ######################################
###################################################################################


output_dump = [];

if __name__=='__main__':
    output_dump.append(_main_test());

