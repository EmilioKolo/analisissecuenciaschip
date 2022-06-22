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
    def __init__(self, genome_name, genome_element='', path_fasta=''):
        # M_seq almacena todos los rangos de secuencias
        # Rangos almacenados en formato chr_n, pos_ini, pos_end
        # Secuencia almacenada en archivos en path a definir
        # Genome element sirve para acelerar cosas pero es necesario si no se corre cargar_promotores_genoma()

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
    # X _consulta_secuencia_fasta(chr_n, pos_ini, pos_end): Devuelve la secuencia consultando en los archivos .fasta
    # - Con funcion cargar_rango() funcional, hacer cargar_bed() y cargar_promotores_genoma(rango)
        # - cargar_bed(archivo): Carga todos los rangos en un archivo de output de ChIP-seq
        # - cargar_promotores_genoma(rango): Carga todos los rangos alrededor de promotores de genes



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
        dir_arch, nom_arch = self._dir_arch(chr_n);
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


    def _consulta_secuencia_fasta(self, chr_n, pos_ini, pos_end):
        # Consulta la secuencia para chr_n entre pos_ini y pos_end de los .fasta

        # Inicializo la variable que se devuelve
        seq = '';

        # Defino el nombre y direccion del archivo a buscar
        dir_arch, nom_arch = self._dir_arch(chr_n);
        # Intento abrir el archivo correspondiente
        try:
            # Abro el archivo para que tire error si no existe
            F = open(dir_arch, 'r');
            F.close();
            # Abro el archivo con seqIO para extraer secuencia
            seq_fasta = self._leer_fasta_chr(dir_arch);

            # Selecciono el rango y lo agrego a seq
            seq = seq_fasta[pos_ini-1:pos_end].seq;
        except:
            logging.warning('Archivo ' + nom_arch + ' no encontrado en la carpeta dada.');
            # Si chr_n no esta presente, defino chr_id en base a chr_n para descargarlo
            chr_id = self._buscar_chrid(chr_n);
            # Si no encuentro chr_id, tiro error
            if chr_id == '':
                logging.error('Cromosoma ' + chr_n + ' no encontrado en lista de IDs para buscar en nucleotide.');
                encontrado = False;
            # Si encuentro chr_id, trato de bajarlo con self._download_chr()
            else:
                encontrado = self._download_chr(chr_n);
            # Si se pudo descargar el archivo, vuelvo a correr _consulta_secuencia_fasta()
            if encontrado:
                seq = self._consulta_secuencia_fasta(chr_n, pos_ini, pos_end);
        return seq


    def _consulta_secuencia_entrez(self, chr_n, pos_ini, pos_end):
        # Confirma la secuencia para chr_n entre pos_ini y pos_end usando Entrez

        # Inicializo la variable que se devuelve
        seq = '';
        # Defino chr_id en base a chr_n
        chr_id = self._buscar_chrid(chr_n);
        # Defino un tiempo para el retry y por si quiero mandar muchas consultas seguidas
        sleep_time = 0.2;
        retry_time = 10;
        # Inicializo un timer para no mandar demasiadas consultas
        time.sleep(sleep_time);
        # Uso try para intentar obtener la secuencia
        try:
            # Obtengo la secuencia directamente con Entrez.efetch()
            handle = Entrez.efetch(db='nucleotide', id=chr_id, rettype='fasta', seq_start=pos_ini, seq_stop=pos_end);
            record = SeqIO.read(handle, 'fasta');
            handle.close();
            seq = record.seq;
        except:
            logging.warning('Problemas con ' + str(chr_n) + ' entre posiciones ' + 
                            str(pos_ini) + ' y ' + str(pos_end) + '.');
            # Si falla, reintento despues de retry_time
            time.sleep(retry_time);
            try:
                # Obtengo la secuencia directamente con Entrez.efetch()
                handle = Entrez.efetch(db='nucleotide', id=chr_id, rettype='fasta', seq_start=pos_ini, seq_stop=pos_end);
                record = SeqIO.read(handle, 'fasta');
                handle.close();
                seq = record.seq;
            except:
                logging.error('Fallo el reintento. Devuelvo secuencia vacia.');
        return seq


    def _leer_fasta_chr(self, dir_arch):
        # Lee el archivo .fasta en dir_arch y devuelve la secuencia
        # Pensado para una secuencia, si hay mas de una tira error y devuelve la primera

        # Abro el archivo con SeqIO
        arch_fasta = SeqIO.parse(dir_arch, 'fasta');
        # Inicializo la variable que se devuelve
        L_out = [];
        # Contador para revisar si hay mas de un registro
        i = 0;
        # Reviso si tiene mas de una secuencia
        for record in arch_fasta:
            i = i + 1;
            L_out.append(record);
        if i != 1:
            logging.warning('Mas de un registro en el .fasta, se devuelve el primero.')
        L_out = L_out[0];
        return L_out


    def _dir_arch(self, chr_n):
        # Devuelve la direccion del archivo correspondiente a chr_n y el nombre del archivo

        # Defino el nombre del archivo a buscar
        nom_arch = self.genome_name + '_' + chr_n + '.fasta';
        if self.path_fasta == '':
            dir_arch = '.\\' + nom_arch;
        else:
            dir_arch = os.path.join(self.path_fasta, nom_arch);
        return dir_arch, nom_arch


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

    # Inicializo la variable que se devuelve
    L_out = [];
    
    # Paso datos de usuario a Entrez
    Entrez.email = 'ekolomenski@gmail.com';
    Entrez.api_key = '9bccbd4949f226a5930868df36c211e8b408';

    # Defino la direccion del .fasta
    path_usado = 'D:\\Archivos doctorado\\Genomas\\';
    # Pruebo inicializar seq_data
    print('>Inicializando base_test.')
    base_test = seq_data('mm9', path_fasta=path_usado); # D:\\Archivos doctorado\\Genomas\\ 
    print('>base_test inicializado.')

    #print('>base_test inicializado. Probando _consulta_secuencia_fasta().')
    #pos_ini = 10000000;
    #pos_end = pos_ini + 100;
    #print('>Secuencia con consulta fasta')
    #print(base_test._consulta_secuencia_fasta('chr1', pos_ini, pos_end));
    #print('>Secuencia con consulta Entrez en clase')
    #print(base_test._consulta_secuencia_entrez('chr1', pos_ini, pos_end));
    #print('>Secuencia con funcion ConsultaSecuencia') ### FUNCION BORRADA DEL ARCHIVO
    #print(ConsultaSecuencia(base_test._buscar_chrid('chr1'), pos_ini, pos_end));

    #print('>base_test inicializado. Cargando rango en chr1.')
    #base_test.cargar_rango('chr1',1000100,1000200);
    #print('>Primer rango cargado en chr1. Cargando segundo rango.')
    #base_test.cargar_rango('chr1',1200100,1200200);
    #print('>Segundo rango cargado en chr1. Cargando rango en chr2.')
    #base_test.cargar_rango('chr2',1000100,1000200);
    #print('>Rango cargado en chr2. Devolviendo dict_range en base_test.')
    #print(base_test.dict_range)


    L_out = base_test;
    return L_out



###################################################################################
######################################## OLD ######################################
###################################################################################


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

