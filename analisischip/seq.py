# Generales
import os
import time
import logging
# Analisis de secuencias y genomas
from Bio import Entrez, SeqIO
from pyensembl import EnsemblRelease


###################################### CLASES #####################################


class seq_data(object):
    '''
Clase para cargar datos de secuencia en un genoma dado
Almacena secuencias en formato chr_n, pos_ini, pos_end
Descarga archivos .fasta con los cromosomas necesarios para las secuencias usadas

FUNCIONES IMPORTANTES:
cargar_rango(chr_n, pos_ini, pos_end): Carga un rango de secuencias en self.dict_range (y hace muchos chequeos)
_chr_check(chr_n): Al cargar, ver si el cromosoma esta presente en self.dict_range o en carpeta path_fasta
_chr_file_check(chr_n): Busca el archivo en carpeta path_fasta y usa _download_chr(chr_n) si no esta
_download_chr(chr_n, retries): Ver consulta_secuencia_chr() en 14-PruebaDescargarChr.py
_buscar_chrid(chr_n): Consigue chrID en base a chr_n, usando self.dict_chrid cargado en __init__
_cargar_dict_chrid(genome_name): Carga self.dict_chrid para pasar de chr_n a chr_id
_agregar_chrid(chr_id, chr_n): Revisa que chr_n no este en self.dict_chrid y agrega chr_id si no esta
_consulta_entrez_chr(chr_id): Consigue el elemento correspondiente al cromosoma con SeqIO
_consulta_secuencia_fasta(chr_n, pos_ini, pos_end): Devuelve la secuencia consultando en los archivos .fasta
cargar_bed(archivo): Carga todos los rangos en un archivo de output de ChIP-seq
cargar_promotores(rango): Carga todos los rangos alrededor de promotores de genes
buscar_sitios_union(L_sitios): Crea y devuelve un elemento seq_data con las posiciones de todos los sitios de union en self.dict_rangos
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


    def _buscar_SU_en_seq(self, seq_union, seq_referencia, pos_ini_ref=0):
        # Busca todas las ocurrencias de seq_union en seq_referencia
        # Devuelve una lista de sitios de union en formato [pos_ini, pos_end, forward]
        # En vez de recorrer dos veces o registrar distintos numeros para forward/reverse, registra pos_ini y pos_end con booleano forward

        # Inicializo la lista de sitios de union
        L_SU = [];
        # Recorro seq_referencia
        for i in range(len(seq_referencia)-len(seq_union)+1):
            # Inicializo curr_SU
            curr_SU = [];
            # Defino la secuencia que se revisa
            curr_seq = str(seq_referencia[i:i+len(seq_union)]).upper();
            # Defino el reverso
            curr_seq_rev = str(self.complemento_secuencia(curr_seq)).upper();
            # Si seq_union es igual a curr_seq, registro pos_ini, pos_end, True
            if str(seq_union).upper() == str(curr_seq).upper():
                curr_SU = [pos_ini_ref+i+1, pos_ini_ref+i+len(seq_union), True];
            # Si seq_union es igual a curr_seq_rev, registro pos_ini, pos_end, False
            elif str(seq_union).upper() == str(curr_seq_rev).upper():
                curr_SU = [pos_ini_ref+i+1, pos_ini_ref+i+len(seq_union), False];
            # Reviso si curr_SU esta registrado
            if len(curr_SU) > 0:
                # Si se encontro un sitio de union, tiene largo mayor a 0 y lo registro en L_SU
                L_SU.append(curr_SU[:]);
        return L_SU


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


    def _complemento(self, N, adn=True):
        # Devuelve el complemento del nucleotido N, para ADN o ARN
    
        # Defino los diccionarios de traduccion de ADN y ARN
        dict_adn = {'T':'A','U':'A','A':'T','C':'G','G':'C','Y':'R','R':'Y','K':'M','M':'K',
                    'W':'W','S':'S','D':'H','H':'D','V':'B','B':'V','N':'N'};
        dict_arn = {'T':'A','U':'A','A':'U','C':'G','G':'C','Y':'R','R':'Y','K':'M','M':'K',
                    'W':'W','S':'S','D':'H','H':'D','V':'B','B':'V','N':'N'};

        # Reviso que el nucleotido a traducir este en las keys de los diccionarios (ambos tienen mismas keys)
        if not (N in dict_adn.keys()):
            # Si N no esta entre las keys, devuelvo nucleotido 'N'
            logging.warning('Nucleotido "' + str(N) + '" no interpretado. Se devuelve N.');
            ret = 'N';
        # Si adn es True, uso dict_adn para la traduccion
        elif adn:
            ret = dict_adn[N];
        # Si adn es False, uso dict_arn para la traduccion
        else:
            ret = dict_arn[N];
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


    def _dir_arch(self, chr_n):
        # Devuelve la direccion del archivo correspondiente a chr_n y el nombre del archivo

        # Defino el nombre del archivo a buscar
        nom_arch = self.genome_name + '_' + chr_n + '.fasta';
        if self.path_fasta == '':
            dir_arch = '.\\' + nom_arch;
        else:
            dir_arch = os.path.join(self.path_fasta, nom_arch);
        return dir_arch, nom_arch


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


    def _obtener_chr(self, contig):
        # Recibe contig de un elemento gen de Entrez y define chr_n en base a eso

        # Lista de contigs que puedo procesar
        lista_contigs = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 
                         '11', '12', '13', '14', '15', '16', '17', '18', '19', 
                         '20', '21', '22', '23', 'X', 'Y', 'M', 'MT'];
        # Reviso que contig este en lista_contigs
        if str(contig).upper() in lista_contigs:
            # Los elementos de la lista de contigs funcionan agregando chr antes del contig
            chr_n = 'chr' + str(contig).upper();
        else:
            logging.warning('Contig "' + str(contig) + '" no se pudo procesar.')
            chr_n = '';
        return chr_n


    def _obtener_genoma(self, genome_version='', organism=''):
        # Asigna un elemento genoma de Ensembl a self.genome
        # Si self.genome_name no es hg19, mm9 ni hg38, intenta buscar usando EnsemblRelease(genome_version, species=organism)

        # Defino el genoma default (por si lo quiero cambiar rapido)
        default_genome = EnsemblRelease(102, species='human');
        # Primero reviso self.genome_name por si es hg19, mm9 o hg38
        if self.genome_name.lower() in ['hg19', 'mm9', 'hg38']:
            # CUIDADO: En init se asume hg19 si no se da uno de estos genome_name
            # Para seguir hay que cambiar init o cambiar manualmente self.genome_name
            if self.genome_name.lower() == 'hg19':
                self.genome = EnsemblRelease(102, species='human');
            elif self.genome_name.lower() == 'mm9':
                self.genome = EnsemblRelease(67, species='mouse');
            elif self.genome_name.lower() == 'hg38':
                self.genome = EnsemblRelease(102, species='mouse');
            else:
                # Esto no deberia pasar
                logging.error('self.genome_name no esta en la lista seleccionada. Se utiliza el genoma por defecto.');
                self.genome = default_genome;
        # Si self.genome no esta cargado y self.genome_name no es hg19, mm9 ni hg38, veo si genome_version o organism tienen info
        elif genome_version != '' or organism != '':
            # Si hay genome_version o organism cargados, veo si ambos estan cargados
            if genome_version != '' and organism != '':
                # Pruebo correr y si falla uso default_genome
                try:
                    self.genome = EnsemblRelease(int(genome_version), species=organism);
                except:
                    logging.error('Falla busqueda de genoma species=' + str(organism)  + '; version=' + str(genome_version) + '. Se utiliza el genoma por defecto.')
                    self.genome = default_genome;
            # Si ninguno esta cargado, tiro error y uso default_genome
            else:
                logging.error('Completar genome_version y organism para buscar en EnsemblRelease. Se utiliza el genoma por defecto.');
                self.genome = default_genome;
        else:
            # Si no hay ningun dato, directamente se asume default_genome
            logging.warning('Sin informacion de genoma. Se utiliza el genoma por defecto.');
            self.genome = default_genome;
        return self


    def _obtener_pos0_forward(self, pos_ini, pos_end, strand, strand_forward='+', strand_reverse='-'):
        # Recibe pos_ini, pos_end y strand de un elemento gen de Entrez
        # En base a eso, define si el gen es forward o reverse y devuelve la posicion del +1
        # Permite cambiar que se define como forward y reverse

        # Reviso strand
        if strand == strand_forward:
            forward = True;
        elif strand == strand_reverse:
            forward = False;
        else:
            logging.error('Strand "' + str(strand) + '" no definida como forward ni reverse. Se usa forward como default.')
            forward = True;
        if forward:
            ret = pos_ini;
        else:
            ret = pos_end;
        return ret, forward


    def buscar_sitios_union_lista(self, L_sitios):
        # Crea y devuelve un elemento seq_data con las posiciones de todos los sitios de union en self.dict_rangos
        # Busca una lista de sitios de union posibles

        # Inicializo el elemento seq_data que se devuelve con los mismos valores de init que el que contiene los rangos
        seq_out = seq_data(self.genome_name, genome_element=self.genome, path_fasta=self.path_fasta);
        # Recorro cada uno de los cromosomas en self.dict_rangos
        for key in self.dict_range.keys():
            L_rangos = self.dict_range[key];
            chr_n = key;
            # Recorro cada rango en L_rangos para chr_n
            for curr_rango in L_rangos:
                # Defino pos_ini y pos_end para curr_rango
                # No reviso la orientacion del rango porque devuelvo el sitio con orientacion incluida
                curr_pos_ini = min(curr_rango[0], curr_rango[1]);
                curr_pos_end = max(curr_rango[0], curr_rango[1]);
                # Busco la secuencia dada en el rango
                seq_rango = self._consulta_secuencia_fasta(chr_n, curr_pos_ini, curr_pos_end);
                # Recorro cada sitio buscado en L_sitios
                for sitio_buscado in L_sitios:
                    # Obtengo una lista de posiciones para los sitios de union encontrados
                    L_SU = self._buscar_SU_en_seq(sitio_buscado, seq_rango, pos_ini_ref=curr_pos_ini-1);
                    # Por cada sitio encontrado, cargo un rango en seq_out
                    for sitio_encontrado in L_SU:
                        # Depende de como devuelvo el sitio en _buscar_SU_en_seq()
                        seq_out.cargar_rango(chr_n, sitio_encontrado[0], sitio_encontrado[1], forward=sitio_encontrado[2]);
        return seq_out


    def cargar_rango(self, chr_n, pos_ini, pos_end, forward=True):
        # Carga un rango de secuencias en self.dict_range (y hace muchos chequeos)
        # No revisa si hay superposicion
        # Forward registra para histogramas en rangos de promotores

        # Uso self._chr_check(chr_n) para revisar que el cromosoma este presente
        chr_presente = self._chr_check(str(chr_n));
        # Si esta presente, trato de agregar pos_ini y pos_end
        if chr_presente:
            # Es posible que sea necesario revisar superposicion, pero inicialmente no lo hago
            # Creo el rango revisando que sean int y pos_ini sea menor que pos_end
            rango_cargado = (min(int(pos_ini), int(pos_end)), max(int(pos_ini), int(pos_end)), forward);
            self.dict_range[chr_n].append(rango_cargado);
        # Si no pude agregar chr_n, tiro error
        else:
            logging.error('Rango no agregado por no poder parsear "' + str(chr_n) + '"')
        return self


    def cargar_promotores(self, rango_promotor, genome_version='', organism=''):
        # Carga todos los rangos del genoma a self.dict_range
        # Registra entre rango_promotor[0] y rango_promotor[1] a partir del +1 de cada gen
        # Usa cargar_rango() con los rangos obtenidos y chequea que el genoma este cargado con _obtener_genoma()

        # Chequeo si self.genome esta cargado
        if self.genome == '':
            # Si no esta cargado, uso self._obtener_genoma() para cargarlo
            self._obtener_genoma(genome_version, organism);

        # Recorro cada gen en el genoma
        for gene_object in self.genome.genes():
            # Solo agarro protein_coding
            if gene_object.biotype == 'protein_coding':
                # Defino pos0 del gen y booleano forward en base a strand y pos_ini/pos_end
                pos0, forward = self._obtener_pos0_forward(gene_object.start, gene_object.end, gene_object.strand);
                # Defino chr_n en base a contig
                chr_n = self._obtener_chr(gene_object.contig);
                # Solo cargo el rango si chr_n se pudo procesar (si no se puede, devuelve string vacio)
                if len(chr_n) > 0:
                    # Si el gen es forward, registro el rango normalmente
                    if forward:
                        self.cargar_rango(chr_n, pos0+rango_promotor[0], pos0+rango_promotor[1],forward=forward);
                    # Si es reverse, tengo que usar rango_promotor al reves y restarlo
                    else:
                        self.cargar_rango(chr_n, pos0-rango_promotor[1], pos0-rango_promotor[0],forward=forward);
        return self


    def complemento_secuencia(self, seq, adn=True):
        # Devuelve el complemento de una secuencia de adn o arn
        # La secuencia del complemento se devuelve al reves que la referencia

        # Inicializo la variable que se devuelve
        ret_seq = '';
        # Recorro seq de atras para adelante
        for i in range(len(seq)):
            # Obtengo el complemento de cada posicion de atras para adelante y las agrego a ret_seq
            ret_seq = ret_seq + self._complemento(seq[-i-1],adn=adn);
        return ret_seq


    def leer_bed(self, nom_bed, path_bed='.\\', col_chr=0, col_ini=1, col_end=2, sep='\t', ext='.bed'):
        # Carga todos los rangos en un archivo .bed con los outputs de ChIP-seq a self.dict_range
        # Usa cargar_rango() con los rangos obtenidos
        # col_chr, col_ini y col_end definen la columna de la matriz donde se encuentran los datos necesarios

        # Defino la direccion del archivo en base a path_bed y nom_bed
        # Si nom_bed no termina en .bed, lo agrego
        if nom_bed[-len(ext):] != ext:
            dir_arch = os.path.join(path_bed, nom_bed + str(ext));
        else:
            dir_arch = os.path.join(path_bed, nom_bed);

        # Intento abrir el archivo .bed en dir_arch
        try:
            # Pruebo abrir el archivo (se pueden agregar mas pruebas)
            F = open(dir_arch, 'r');
            F.close();
            # Registro que el archivo existe
            arch_existe = True;
        except:
            # Si tira error, registro que el archivo no existe
            arch_existe = False;
        # Si el archivo existe, lo abro propiamente
        if arch_existe:
            with open(dir_arch, 'r') as F:
                # Recorro cada linea del archivo
                for curr_line in F:
                    # Hago una lista de la linea actual
                    L = curr_line.rstrip().split(sep);
                    # Reviso que L tenga largo para los valores de busqueda de columnas
                    if len(L) > max(col_chr, col_ini, col_end):
                        # Defino chr, pos_ini y pos_end en base a los valores de columna dados
                        chr_bed = L[col_chr];
                        pos_ini = int(L[col_ini]);
                        pos_end = int(L[col_end]);
                        # Si chr_bed empieza con chr, lo uso como chr_n
                        if chr_bed[:3] == 'chr':
                            self.cargar_rango(chr_bed, pos_ini, pos_end);
                        # Si chr_bed es un string corto, reviso con _obtener_chr()
                        elif len(chr_bed) <= 2:
                            chr_n = self._obtener_chr(chr_bed);
                            # Si _obtener_chr() devuelve un string no vacio, uso cargar_rango()
                            if len(chr_n) > 0:
                                self.cargar_rango(chr_n, pos_ini, pos_end);
                            # Si devuelve un string vacio, tiro error
                            else:
                                logging.error('No se pudo interpretar "' + str(chr_bed) + '".');
                        # Si chr_bed es un string largo que no empieza con 'chr', tiro error
                        else:
                            logging.error('No se pudo interpretar "' + str(chr_bed) + '".');
                    # Si L no tiene el largo, tiro error
                    else:
                        logging.warning('Fila demasiado corta. Lista obtenida: ' + str(L) + '.')
        # Si el archivo no se pudo abrir, tiro error
        else:
            logging.error('No se pudo abrir el archivo "' + str(nom_bed) + '" en "' + str(path_bed) + '".');
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
    print('>base_test inicializado. ')


    #print('>base_test inicializado. Cargando rangos para probar busqueda de sitios de union.')
    #base_test.cargar_rango('chr1',10000100,10000200);
    #base_test.cargar_rango('chr1',10200100,10200200);
    #base_test.cargar_rango('chr1',10001000,10002000);
    #print('>Rangos cargados. Buscando sitios de union.')
    #for i in base_test.dict_range['chr1']:
    #    print(i)
    #    print(base_test._consulta_secuencia_fasta('chr1',i[0],i[1]))
    #L_sitios = ['AAGTG'];
    #sitios_test = base_test.buscar_sitios_union_lista(L_sitios);
    #print('>Rangos de sitios de union encontrados.')
    #print(sitios_test.dict_range)
    #for i in sitios_test.dict_range['chr1']:
    #    print(i)
    #    print(sitios_test._consulta_secuencia_fasta('chr1',i[0],i[1]))

    #print('>base_test inicializado. Probando busqueda de seq en seq.')
    #seq_ref = 'ATATTACGATCGT';
    #seq_busq = 'TCGT';
    #L_SU = base_test._buscar_SU_en_seq(seq_busq, seq_ref);
    #print('Secuencia de referencia:')
    #print(seq_ref)
    #print('Secuencia buscada:')
    #print(seq_busq)
    #print('Posiciones encontradas:')
    #for SU in L_SU:
    #    print(SU)

    #print('>base_test inicializado. Probando carga de archivo .bed.')
    #nom_bed = 'Dupays2015';
    #path_bed = 'D:\\Dropbox\\Doctorado\\3-Genes transactivados rio abajo\\0-Fuentes\\Papers ChIP-seq\\';
    #print('>Probando leer_bed con "' + nom_bed + '" en "' + path_bed + '".')
    #base_test.leer_bed(nom_bed, path_bed);
    #print('>Archivo leido y cargado. Mostrando dict_range.')
    #for key in base_test.dict_range.keys():
    #    print(key)
    #    print(base_test.dict_range[key])

    #print('>base_test inicializado. Inicializando revision de genoma.')
    #base_test.cargar_promotores([-1500, 1500]);
    #print('>Carga de promotores finalizada. Mostrando dict_range.')
    #for key in base_test.dict_range.keys():
    #    print(key)
    #    print(base_test.dict_range[key])

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

    L_out = base_test.dict_range;
    return L_out



###################################################################################
######################################## OLD ######################################
###################################################################################


def buscar_en_secuencia(busq, seq):
    # Busca una secuencia "busq" en una secuencia mas larga "seq"
    # Devuelve la posicion en la que se encontro busq y True si se encuentra busq
    # Devuelve len(seq)-len(busq)+1 y False si no se encuentra busq
    pos = 0;
    encontrado = False;
    # Recorro una por una las posiciones de seq hasta encontrar busq
    while pos < (len(seq)-len(busq)+1) and (not encontrado):
        if seq[pos:pos+len(busq)] == busq:
            encontrado = True;
        else:
            pos += 1;
    return pos, encontrado


def buscar_en_seq_ambas_direcciones(busq, seq):
    # Funcion que usa buscar_en_secuencia() para buscar todas las ocurrencias de busq en seq
    # Busca en ambas direcciones y registra varias ocurrencias
    L_pos = [];
    loop_1 = True;
    curr_pos = 0;
    seq_pos = 0;
    # Primero corro buscar_en_secuencia() para busq y seq normalmente
    while loop_1:
        curr_pos, loop_1 = buscar_en_secuencia(busq, seq[seq_pos:]);
        if loop_1:
            L_pos.append(seq_pos+curr_pos);
            seq_pos = seq_pos + curr_pos + 1;
    # Despues vuelvo a correr con la secuencia reversa
    rev_seq = seq_data.complemento_secuencia(seq, adn=True);
    loop_1 = True;
    curr_pos = 0;
    seq_pos = 0;
    while loop_1:
        curr_pos, loop_1 = buscar_en_secuencia(busq, rev_seq[seq_pos:]);
        if loop_1:
            L_pos.append(-(seq_pos+curr_pos));
            seq_pos = seq_pos + curr_pos + 1;
    return L_pos


def buscar_en_seq_2dir_unificado(busq, seq):
    # Funcion que usa buscar_en_secuencia() para buscar todas las ocurrencias de busq en seq
    # Busca en ambas direcciones y registra varias ocurrencias

    L_pos = [];
    loop_1 = True;
    curr_pos = 0;
    seq_pos = 0;
    # Primero corro buscar_en_secuencia() para busq y seq normalmente
    while loop_1:
        curr_pos, loop_1 = buscar_en_secuencia(busq, seq[seq_pos:]);
        if loop_1:
            n = seq_pos+curr_pos;
            L_pos.append(n-len(seq));
            seq_pos = seq_pos + curr_pos + 1;
    # Despues vuelvo a correr con la secuencia reversa
    rev_seq = seq_data.complemento_secuencia(seq, adn=True);
    loop_1 = True;
    curr_pos = 0;
    seq_pos = 0;
    while loop_1:
        curr_pos, loop_1 = buscar_en_secuencia(busq, rev_seq[seq_pos:]);
        if loop_1:
            n = seq_pos+curr_pos;
            L_pos.append(-n-len(busq));
            seq_pos = seq_pos + curr_pos + 1;
    return L_pos


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

