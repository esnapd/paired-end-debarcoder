# Handle Barcode Information
#
#
#
#
#
import csv

from cytoolz.dicttoolz import assoc
from cytoolz.functoolz import thread_first

from schema import Schema, And, Or

###################################################################
## Schemas
barcodeSchema = Schema(
    {"barcode":         And(str, lambda s:  12 <= len(s) <= 30),
     "forward_barcode": And(str, lambda s:  6 <= len(s) <= 15),
     "forward_spacer":  And(str, lambda s:  0 <= len(s) <= 10),
     "forward_primer":  And(str, lambda s:  12 <= len(s) <= 35),
     "reverse_barcode": And(str, lambda s:  6 <= len(s) <= 15),
     "reverse_spacer":  And(str, lambda s:  0 <= len(s) <=10),
     "reverse_primer":  And(str, lambda s:  12 <= len(s) <= 35)
     })

fastadataSchema = Schema(
    {"forward_id":       str,
     "forward_desc":     str,
     "forward_sequence": str,
     "reverse_id":       str,
     "reverse_desc":     str,
     "reverse_sequence": str,
     "sample" :          Or(str, None),
     "barcode":          And(str, lambda s:  12 < len(s) < 30),
     "barcode_distance": int,
     "tooshort":         bool,
     "spacermismatch":   bool,
     })


def process_barcodefile(file, barcodelength, checkbarcodes=True):
    "Take a barcode file and return the barcode"

    data = {}
    lines = open(file,'r')
    reader = csv.reader(lines.readlines(), delimiter='\t')

    for idx, line in enumerate(reader):
        if idx > 0 and line:   # don't use header and don't use empty lines
            if len(line) == 1:
                raise ValueError("There is a problem with your barcode file. Is it tab-delimited?")
            try:
                if len(line) == 8:
                    sample, barcode, forward_barcode, forward_spacer, forward_primer, \
                    reverse_barcode, reverse_spacer, reverse_primer  = line
                elif len(line) > 8:
                    sample, barcode, forward_barcode, forward_spacer, forward_primer, \
                    reverse_barcode, reverse_spacer, reverse_primer, *othercols = line
                else:
                    raise ValueError("Barcode File must have a minimum of 8 data columns")
            except:
                raise ValueError("Barcode File must have a minimum of 8 data columns")

            # validate the barcode data
            barcodedata =  {"barcode":         barcode,
                            "forward_barcode": forward_barcode,
                            "forward_spacer":  forward_spacer,
                            "forward_primer":  forward_primer,
                            "reverse_barcode": reverse_barcode,
                            "reverse_spacer":  reverse_spacer,
                            "reverse_primer":  reverse_primer}

            if checkbarcodes:
                barcodedata = barcodeSchema.validate(barcodedata)

            data[sample] = barcodedata

    #check data
    if checkbarcodes:
        assert(data != {})
        for k,v in data.items():
            # check barcode lengths
            if not len(v['barcode']) == barcodelength:
                raise ValueError("Barcode {}, of sample {} is not of expected length {}".format(v['barcode'], k, barcodelength))
            #check forward and reverse barcodes
            assert(v['forward_barcode'] + v['reverse_barcode'] == v['barcode'])

    barcodes = [v['barcode'] for k,v in data.items()]
    if not len(barcodes) == len(set(barcodes)):
        raise ValueError("Barcode Values are not Unique. Please Check your Barcoding File")

    return data

def check_barcode(fastadict, barcodedict, barcodelength, maxdistance):
    "check for barcode and update sample data"

    samplematch      = None
    barcodedata      = None
    spacermismatch   = False
    barcode_distance = 0
    halfbarcode      = int(barcodelength/2)
    fseq    = fastadict['forward_sequence']
    rseq    = fastadict['reverse_sequence']
    barcode = fseq[:halfbarcode] + rseq[:halfbarcode]

    #check for perfect match first:
    for	sample, samplebarcodedict in barcodedict.items():
        if samplebarcodedict['barcode'] == barcode:
            samplematch = sample
            barcodedata = samplebarcodedict

    #if not choose closest
    if not samplematch:
        for	sample, samplebarcodedict in barcodedict.items():
            hdist = hamdist(samplebarcodedict['barcode'], barcode)
            if hdist <= maxdistance:
                barcode_distance = hdist
                samplematch = sample
                barcodedata = samplebarcodedict

    # trim the sequences after checking the spacer sequence between the barcode and the primer
    fseq = fseq[halfbarcode:]
    rseq = rseq[halfbarcode:]

    if barcodedata is not None:
        forward_spacer = barcodedata['forward_spacer']
        reverse_spacer = barcodedata['reverse_spacer']

        if fseq.startswith(forward_spacer):
            fseq = fseq[len(forward_spacer):]
        else:
            fseq = fseq[len(forward_spacer):]
            spacermismatch = True
        if rseq.startswith(reverse_spacer):
            rseq = rseq[len(reverse_spacer):]
        else:
            rseq = rseq[len(reverse_spacer):]
            spacermismatch = True

    # return updated values
    return thread_first(fastadict,
                        (assoc, "sample", samplematch),
                        (assoc, "spacermismatch", spacermismatch),
                        (assoc, "barcode", barcode),
                        (assoc, "barcode_distance", barcode_distance),
                        (assoc, "forward_sequence", fseq),
                        (assoc, "reverse_sequence", rseq))


def hamdist(str1, str2):
   "Count the # of differences between equal length strings str1 and str2"
   diffs = 0
   for ch1, ch2 in zip(str1, str2):
       if ch1 != ch2:
           diffs += 1
   return diffs