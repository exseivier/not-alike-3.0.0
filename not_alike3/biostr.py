#!/usr/bin/env python3
import os
import sys
import ctypes as ct


#
#   BIOSTRINGS PYTHON MODULE
#   
#   AUTHOR: JAVIER MONTALVO-ARREDONDO
#   CONTACT: buitrejma@gmail.com
#   UNIVERSIDAD AUTONOMA AGRARIA ANTONIO NARRO
#   DEPARTAMENTO DE CIENCIAS BASICAS.
#   BASIC SCIENCES DEPARTMENT.
#

#  TODO:    Problem is the version of python in this path string.
#           Has to be changed to correct the version depending on the current version
#           the final user is using.
#           Suggestion
#               PyVersion = "python" + ".".join(sys.version.split(" ")[0].split(".")[:2])
#               pathString = f'lib/{PyVersion}/site-packages/not_alike3/biostruct'
#               file_path = os.path.join(sys.prefix, pathString)
#   XXX: Have to test it. No tested yet.

PyVersion = "python" + ".".join(sys.version.split(" ")[0].split(".")[:2])
pathString = f'lib/{PyVersion}/site-packages/not_alike3/biostruct'
file_path = os.path.join(sys.prefix, pathString)

libc = ct.CDLL(f'{file_path}/libdnah.so')

class DNA(ct.Structure):
    """
        DNA structure which holds the sequence,
        IDs, start and end position where it starts and ends.
    """

    _fields_ = [
                ('seq', ct.c_char_p),
                ('ids', ct.POINTER(ct.c_char_p)),
                ('start', ct.POINTER(ct.c_int)),
                ('end', ct.POINTER(ct.c_int)),
                ('seq_len', ct.POINTER(ct.c_int)),
                ('hide', ct.POINTER(ct.c_int))
                ]

class lkdList(ct.Structure):
    """
        Linked list structure to hold headers.
    """
    pass

lkdList._fields_ = [
                ('header', ct.c_char_p),
                ('next', ct.POINTER(lkdList))
                ]

libc.loadDNASeqs.restype = ct.POINTER(DNA)
libc.loadDNASeqs.argtypes = {ct.c_char_p}

libc.freeLkdList.argtypes = [ct.POINTER(lkdList)]

libc.freeDNA.argtypes = [ct.POINTER(DNA)]

libc.writeNoHideToFile.argtypes =[ct.POINTER(DNA)]

libc.splitBioString.restype = ct.POINTER(DNA)
libc.splitBioString.argtypes = [ct.POINTER(DNA), ct.c_int, ct.c_int]

libc.loadLines_lkdList.restype = ct.POINTER(lkdList)
libc.loadLines_lkdList.argtypes = [ct.c_char_p]

libc.filterBioseq.argtypes = [ct.POINTER(DNA), ct.POINTER(lkdList)]

libc.sampleSeqs.argtypes = [ct.POINTER(DNA), ct.c_int, ct.c_char_p]

class DNAstuffs():

    a_ptr = None
    is_empty = True
    size = None
    step = None
 
    def freeDNA(self):
        """
            Free allocated memory
        """
        if self.is_empty == True:
            print("Object is empty")
            return None
        try:
            libc.freeDNA(self.a_ptr)
            print("Free dynamic allocated memory for DNA!")
        except Exception as e:
            print(f'Something went worng!: {e}')

        self.is_empty = True

    def getSeq(self, ids):
        """
            Retrieves a sequence by id
        """
        if self.is_empty == True:
            print("Object is empty")
            return None
        i = 0
        while self.a_ptr.contents.ids[i] != None:
            seq_id = self.a_ptr.contents.ids[i].decode()
            start = self.a_ptr.contents.start[i]
            end = self.a_ptr.contents.end[i]
            print(f'{start}, {end}')
            seq = self.a_ptr.contents.seq.decode()[start : (end + 1)]
            if seq_id == ids:
                return f'{seq_id}\n{seq}\n'
            i += 1

        return "ID does not exist!"

    def getObject(self):
        if self.is_empty == True:
            print("Object is empty")
            return None
        else:
            return self.a_ptr

    def writeNoHideToFile(self, outfile):
        """
            Writes no-hided sequences to file.
        """
        try:
            libc.writeNoHideToFile(self.a_ptr, outfile.encode())
#            return True
        except Exception as e:
            print(f'Gess what! Something went wrong in writeNoHideToFile!')
            print(f'Because of: {e}')
#            return False

    def filterBioseq(self, lkdList_heads):
        """
            Drop off those sequences whose headers are in the linked list of headers.
            So that means those sequences had at least a hit in BLASTn search.
        """
        try:
            libc.filterBioseq(self.a_ptr, lkdList_heads.getObject())
            return True
        except Exception as e:
            print(f'Guess what! Somthing went wrong in filterBioseq!')
            print(f'Beacuse of: {e}')

    def sampleSeqs(self, perc, outfile = 'sampleout.fasta'):

        assert perc >= 10, f'A number greater or equal to 10 is expected, got: {perc}'
    
        libc.sampleSeqs(self.a_ptr, perc, outfile.encode())


class loadDNASeqs(DNAstuffs):
    """
        loadDNASeqs objcet
    """
    def __init__(self, filename):
        """
            Initializes the object
        """
        self.a_ptr = ct.POINTER(DNA)
        self.a_ptr = libc.loadDNASeqs(filename.encode())
        self.is_empty = False

    def __str__(self):
        """
            Orints info from object
        """
        if self.a_ptr != None:
            return "Object loaded!"
        else:
            return "Object unloaded!"


class splitBioString(DNAstuffs):
    """
        Split BioString class.
    """

    def __init__(self, seqs, size, step):

        self.a_ptr = ct.POINTER(DNA)
        self.a_ptr = libc.splitBioString(seqs.getObject(), size, step)
        self.is_empty = False

    def __str__(self):

        if self.is_empty == True:
            return "Object is empty!"
        else:
            return "Object is loaded!"

class loadLkdList():
    """
        Loads a list of headers and stores them into lkdList object.
    """
    lkdl = None
    is_empty = True

    def __init__(self, filename):
        """
            Initializing
        """
        self.lkdl = ct.POINTER(lkdList)
        self.lkdl = lkdList()
        
        self.lkdl = libc.loadLines_lkdList(filename.encode())
        self.is_empty = False

    def __str__(self):
        """
            Shows on screen if object is loaded or not
        """
        if self.is_empty == False:
            return "Object is loaded!"
        else:
            return "Object is empty!"


    def getObject(self):
        """
            Return lkdList object.
        """
        if self.is_empty == False:
            return self.lkdl
        else:
            print("Object is empty!")
            return None

    def findIt(self, head):
        """
            Finds and returns the line where line == head. 
        """
        tmp_lkdl = ct.POINTER(lkdList)
        tmp_lkdl = self.lkdl
        while True:
            try:
                header = tmp_lkdl.contents.header
                if header.decode() == head:
                    return True
                tmp_lkdl = tmp_lkdl.contents.next
            except Exception as e:
                return False

    def freeLkdList(self):
        """
            Free dynamic allocated memory for a lkdList object instance.
        """
        libc.freeLkdList(self.lkdl)
        print("Free dynamic allocated memory for this lkdList object")


def main():
    """
        Main program
    """

    seqs = loadDNASeqs('examples/seqs.fas')

#    print(seqs.getSeq('>Seq1'))
#    print(seqs.getObject())

    sptSeqs = splitBioString(seqs, 20, 10)

    seqs.freeDNA()

#    print(sptSeqs.getSeq(">Seq1_0"))
#    print(sptSeqs.getSeq(">Seq1_1"))
#    print(seqs.getSeq('>Seq1'))
#    print(seqs.getObject())

    headers = loadLkdList('examples/heads.txt')
    
#    print(headers.findIt('>Seq4'))
    
    sptSeqs.filterBioseq(headers)
    sptSeqs.writeNoHideToFile()
    sptSeqs.sampleSeqs(10)

    headers.freeLkdList()
    sptSeqs.freeDNA()


if __name__ == '__main__':
    main()
