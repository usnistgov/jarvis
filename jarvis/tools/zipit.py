import os, sys, zipfile


def ZipDir(inputDir, outputZip, contents=[]):
    """
    Zip up a directory and preserve symlinks and empty directories
    """
    zipOut = zipfile.ZipFile(outputZip, "w", compression=zipfile.ZIP_DEFLATED)
    tmp = contents
    rootLen = len(os.path.dirname(inputDir))

    def _ArchiveDirectory(parentDirectory):
        # contents = os.listdir(parentDirectory)
        contents = (
            tmp
        )  # ['init.mod', 'potential.mod', 'in.elastic', 'data',  'log.lammps', 'restart.equil','data0' ]
        # store empty directories
        if not contents:
            # http://www.velocityreviews.com/forums/t318840-add-empty-directory-using-zipfile.html
            archiveRoot = parentDirectory[rootLen:].replace("\\", "/").lstrip("/")
            zipInfo = zipfile.ZipInfo(archiveRoot + "/")
            zipOut.writestr(zipInfo, "")
        for item in contents:
            fullPath = os.path.join(parentDirectory, item)
            if os.path.isdir(fullPath) and not os.path.islink(fullPath):
                _ArchiveDirectory(fullPath)
            else:
                archiveRoot = fullPath[rootLen:].replace("\\", "/").lstrip("/")
                if os.path.islink(fullPath):
                    # http://www.mail-archive.com/python-list@python.org/msg34223.html
                    zipInfo = zipfile.ZipInfo(archiveRoot)
                    zipInfo.create_system = 3
                    # long type of hex val of '0xA1ED0000L',
                    # say, symlink attr magic...
                    # zipInfo.external_attr = 0777 << 16L
                    zipInfo.external_attr = "2716663808L"
                    zipOut.writestr(zipInfo, os.readlink(fullPath))
                else:
                    zipOut.write(fullPath, archiveRoot, zipfile.ZIP_DEFLATED)

    _ArchiveDirectory(inputDir)

    zipOut.close()
