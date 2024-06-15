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
#   XXX: Tested OK.

PyVersion = "python" + ".".join(sys.version.split(" ")[0].split(".")[:2])
pathString = f'lib/{PyVersion}/site-packages/not_alike3/biostruct'
file_path = os.path.join(sys.prefix, pathString)

libc = ct.CDLL(f'{file_path}/libdnah.so')
#libc = ct.CDLL(f'/home/jima/projects/Dev/python_apps/biostrc/biostr/biostruct/libdnah.so')

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
                ('hide', ct.POINTER(ct.c_int)),
                ('lhand', ct.c_int),
                ('rhand', ct.c_int),
                ('num_items', ct.c_int),
                ('memAllocSeq', ct.c_int),
                ('numItemsAlloc', ct.c_int),
                ('size', ct.c_int)
                ]

class Node(ct.Structure):
    """
        Node structure
    """
    pass

Node._fields_ = [
                ('header', ct.c_char_p),
                ('next', ct.POINTER(Node))
                ]

class lkdList(ct.Structure):
    """
        Linked list structure to hold headers.
    """
    pass

lkdList._fields_ = [
                    ('head', ct.POINTER(Node)),
                    ('tail', ct.POINTER(Node))
                    ]
# 
#   Biostr library
#

libc.Biostr_init.argtypes = [ct.POINTER(DNA)]
libc.Biostr_init_custom.argtypes = [ct.POINTER(DNA), ct.c_ulonglong]
libc.Biostr_load_file.argtypes = [ct.POINTER(DNA), ct.c_char_p]

libc.Biostr_slice.restype = ct.POINTER(DNA)
libc.Biostr_slice.agrtypes = [ct.POINTER(DNA), ct.c_int, ct.c_int]

libc.Biostr_write_to_file.argtypes = [ct.POINTER(DNA), ct.c_char_p]

libc.Biostr_filter.argtypes = [ct.POINTER(DNA), ct.POINTER(lkdList)]

libc.Biostr_free.argtypes = [ct.POINTER(DNA)]

libc.Biostr_print.argtypes = [ct.POINTER(DNA)]

libc.Biostr_sample.argtypes = [ct.POINTER(DNA), ct.c_int, ct.c_char_p]

#
#   Linked list library
#

libc.create_linked_list.restype = ct.POINTER(lkdList)

libc.create_node.restype = ct.POINTER(Node)
libc.create_node.argtypes = [ct.c_char_p]

libc.add_node.argtypes = [ct.POINTER(lkdList), ct.c_char_p]

libc.free_linked_list.argtypes = [ct.POINTER(lkdList)]

libc.free_node.argtypes = [ct.POINTER(Node)]

libc.insert_at_position.argtypes = [ct.POINTER(lkdList), ct.c_char_p, ct.c_int]

libc.print_linked_list.argtypes = [ct.POINTER(lkdList)]

libc.load_to_linkedList.argtypes = [ct.POINTER(lkdList), ct.c_char_p]

#
#   FileIO libary
#

libc.open_file.restype = ct.c_char_p
libc.open_file.argtypes = [ct.c_char_p]

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
            libc.Biostr_free(self.a_ptr)
            print("Free dynamic allocated memory for DNA!")
            self.is_empty = True
            self.a_ptr = None
        except Exception as e:
            print(f'Something went worng!: {e}')

    def getSeq(self, ids):
        """
            Retrieves a sequence by id
        """
        if self.is_empty:
            print("Object is empty")
            return None
        i = 0
        while self.a_ptr.contents.ids[i] != None:
            seq_id = self.a_ptr.contents.ids[i].decode()
            if seq_id == ids:
                start = self.a_ptr.contents.start[i]
                end = self.a_ptr.contents.end[i]
                seq = self.a_ptr.contents.seq.decode()[start : (end + 1)]
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
            libc.Biostr_write_to_file(self.a_ptr, outfile.encode())
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
            libc.Biostr_filter(self.a_ptr, lkdList_heads.getObject())
            return True
        except Exception as e:
            print(f'Guess what! Somthing went wrong in filterBioseq!')
            print(f'Beacuse of: {e}')

    def sampleSeqs(self, perc, outfile = 'sampleout.fasta'):

        assert perc >= 10, f'A number greater or equal to 10 is expected, got: {perc}'
        
        #   XXX: This function was copied from dnah.c and compiled in new libdnah.so
        #   XXX: Tested before production.
        libc.Biostr_sample(self.a_ptr, perc, outfile.encode())

    def printIt(self):
        
        if self.is_empty:
            print(f'Object is empty')
        else:
            libc.Biostr_print(self.a_ptr)

class loadDNASeqs(DNAstuffs):
    """
        loadDNASeqs object
    """
    def __init__(self, filename):
        """
            Initializes the object
        """
        if os.path.exists(filename) and os.path.isfile(filename):
            self.a_ptr = DNA()
            self.a_ptr = ct.pointer(self.a_ptr)
            libc.Biostr_init(self.a_ptr)
            libc.Biostr_load_file(self.a_ptr, filename.encode())
            self.is_empty = False
        else:
            print(f'FileIOError: Object creation failure: {filename} was not found')

    def __str__(self):
        """
            Prints info from object
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
        
        self.a_ptr = DNA()
        self.a_ptr = ct.pointer(self.a_ptr)
        self.a_ptr = libc.Biostr_slice(seqs.getObject(), size, step)
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

        if os.path.exists(filename) and os.path.isfile(filename):
            self.lkdl = lkdList()
            self.lkdl = ct.pointer(self.lkdl)
            self.lkdl = libc.create_linked_list()
            libc.load_to_linkedList(self.lkdl, filename.encode())
        else:
            print(f"FileIOError: Object creation failure: {filename} was not found")

        self.is_empty = False

    def __str__(self):
        """
            Shows on screen if object is loaded or not
        """
        if self.is_empty == False:
            return "Object is loaded!"
        else:
            return "Object is empty!"

    def printIt(self):

        if self.is_empty:
            print("Object is empty")
        else:
            libc.print_linked_list(self.lkdl)


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
        tmp_lkdl = lkdList()
        tmp_lkdl.head = self.lkdl.contents.head
        tmp_lkdl.tail = self.lkdl.contents.tail
        while True:
            try:
                if tmp_lkdl.head.contents.header.decode() == head:
                    return True
                tmp_lkdl.head = tmp_lkdl.head.contents.next
            except Exception as e:
                print(f"FindError: {e}")
                return False

    def freeLkdList(self):
        """
            Free dynamic allocated memory for a lkdList object instance.
        """
        if self.is_empty:
            print("Linked list object is already empty!")
        else:
            libc.free_linked_list(self.lkdl)
            self.lkdl = None
            self.is_empty = True
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
