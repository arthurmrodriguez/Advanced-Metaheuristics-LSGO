#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Recursively find all *.info and *.pickle files within a directory."""

import os
import warnings
#import zipfile
#import tarfile

# Initialization
def main(directory=os.getcwd(), verbose=True):
    """Lists *.info and *.pickle files recursively in a given directory."""

    filelist = list()

    #~ if directory.endswith('.zip'):
        #~ archive = zipfile.ZipFile(directory)
        #~ for elem in archive.namelist():
            #~ if elem.endswith('.info'):
                #~ (root,elem) = os.path.split(elem)
                #~ filelist = IndexFile(root,elem,archive)
    #~ if directory.find('.tar') != -1:
        #~ archive = tarfile.TarFile(directory)
        #~ for elem in archivefile.namelist():
            #~ if elem.endswith('.info'):
                #~ (root,elem) = os.path.split(elem)
                #~ filelist = IndexFile(root,elem,archive)
    #~ else:

    # Search through the directory directory and all its subfolders.
    for root, dirs, files in os.walk(directory):
        if verbose:
            print 'Searching in %s ...' % root

        for elem in files:
            if elem.endswith('.info') or elem.endswith('.pickle'):
                filelist.append(os.path.join(root, elem))

    if verbose:
        print 'Found %d file(s)!' % (len(filelist))
    if not filelist:
        warnings.warn('Could not find any file of interest in %s!' % root)
    return filelist

if __name__ == '__main__':
    main()
