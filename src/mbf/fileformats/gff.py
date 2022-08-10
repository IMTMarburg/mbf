from mbf_externals.util import to_bytes
from .util import open_file
from .bed import normalize_strand, read_bed


def _mapGFF(row):
    res = {
        'seqname': row[0],
        'source': row[1],
        'feature': row[2],
        'start': int(row[3]) if row[3] != '.' else 0,
        'end': int(row[4]),
        'score': row[5],
        'strand': normalize_strand(row[6]),
        'frame': row[7],
        'attributes': {},
        'comment': ''
    }
    attributes = row[8]
    comment = ''
    if attributes.find(b'#') != -1:
        comment = attributes[attributes.find(b'#') + 1:]
        attributes = attributes[:attributes.find(b'#')]
    res['comment'] = comment
    for x in attributes.split(b';'):
        if x.strip():
            if b'=' in x:
                x = x.split(b"=")
            else:
                x = x.split()
            res['attributes'][x[0].decode('utf-8')] = [y.decode('utf-8') for y in x[1:]]
    return res


def gffToDict(filename, comment_char=None):
    comment_char = to_bytes(comment_char)
    o = open_file(filename)
    rows = o.readlines()
    rows = (x.strip() for x in rows)
    rows = (x.split(b"\t") for x in rows if x and (comment_char is None or x[0] != comment_char[0]))
    res = [_mapGFF(x) for x in rows]
    return res


def dictsToGFF(gffDicts, filename):
    rows = []
    for entry in gffDicts:
        if 'score' in entry:
            score = entry['score']
        else:
            score = '.'
        if 'attributes' in entry:
            attributes = entry['attributes']
            if isinstance(attributes, dict):
                attributes_str = []
                for key, value in attributes.items():
                    attributes_str.append("%s=%s" % (key, value))
                attributes_str = ";".join(attributes_str)
        else:
            attributes_str = ''
        if entry['strand'] not in ('+', '-', '.'):
            raise ValueError("invalid strand: %s" % entry['strand'])
        try:
            entry['frame'] = int(entry['frame'])
            if entry['frame'] not in (0, 1, 2):
                raise ValueError()
        except ValueError:
            raise ValueError("invalid frame: %s " % entry['frame'])
        row = "\t".join((entry['seqname'], entry['source'], entry['feature'], str(int(entry['start'])), str(int(entry['end'])),
                         str(score), entry['strand'], str(entry['frame']), attributes_str))
        rows.append(row)
    o = open(filename, 'wb')
    o.write("\n".join(rows))
    o.close()


class GFF3:

    def escape(self, str):
        if str is None:
            return '.'
        escape = '"\t\n\r=;'
        for k in escape:
            str = str.replace(k, '%' + '%2x' % ord(k))
        return str

    def format_attributes(self, attributes):
        if attributes is None:
            return '.'
        if isinstance(attributes, dict):
            attributes = list(attributes.items())
        valid_attributes = [ 'Name', 'Alias', 'Parent', 'Target', 'Gap', 'Derives_from', 'Note', 'Dbxref']
        res = []
        for id, value in attributes:
        #    if not id in valid_attributes and id != id.lower(): #lower case names are not reserved
         #       raise ValueError("Not a valid tag: %s" % id)
            res.append('%s=%s' % (self.escape(id), self.escape(value)))
        return ";".join(res)


    def dump_row(self, file_handle, seqid = None, source = None, type = None, start = None, end = None, score = None, strand = None, phase = None, attributes = None):
        file_handle.write("\t".join((
            self.escape(seqid),
            self.escape(source),
            self.escape(type),
            self.escape(start),
            self.escape(end),
            self.escape(score),
            self.escape(strand),
            self.escape(phase),
            self.format_attributes(attributes))))
        file_handle.write("\n")


def bed_to_gff(input_filename_or_handle, output_filename_or_handle, feature_type, source = None, show_progress= False):
    entries = read_bed(input_filename_or_handle, report_progress = show_progress)
    output_file_handle = open_file(output_filename_or_handle, 'wb')
    gff = GFF3()
    for e in entries:
        gff.dump_row(output_file_handle, e.refseq, source=source, type=feature_type, start = str(e.position), end = str(e.position + e.length), strand='.')
    output_file_handle.close()
