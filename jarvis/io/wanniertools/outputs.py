"""Class for analyzing  WT.out file."""


class WTOut(object):
    """Construct WT.out related object."""

    def __init__(self, path=""):
        """Initialize with path to WT.out file."""
        self.path = path

    def get_z2_index(self):
        """Get Z2 index."""
        f = open(self.path, "r")
        lines = f.read().splitlines()
        f.close()
        strng = 0
        for ii, i in enumerate(lines):
            if "# z2 number for 6 planes" in i:
                for j in range(6):
                    tmp = lines[ii + j + 1].split()

                    if (
                        tmp[0] == "k1=0.0,"
                        or tmp[0] == "k2=0.0,"
                        or tmp[0] == "k3=0.0,"
                        or tmp[0] == "k3=0.5,"
                        or tmp[0] == "k2=0.5,"
                        or tmp[0] == "k1=0.5,"
                    ):
                        val = tmp[3]
                        strng = strng + float(val)
                    if tmp[0] == "k1=0.5,":
                        weak1 = tmp[3]
                    if tmp[0] == "k2=0.5,":
                        weak2 = tmp[3]
                    if tmp[0] == "k3=0.5,":
                        weak3 = tmp[3]
        index = (
            str(int(strng % 2))
            + str(";")
            + str(weak1)
            + str(",")
            + str(weak2)
            + str(",")
            + str(weak3)
        )
        return index

    def get_chern_number(self):
        """Get Chern index."""
        f = open(self.path, "r")
        lines = f.read().splitlines()
        f.close()
        chrn = []
        try:
            for j, jj in enumerate(lines):
                if "# Chern number for 6 planes" in jj:
                    for k in range(1, 7):
                        tmp = float(lines[j + k].split(":")[1])
                        if tmp not in chrn:
                            chrn.append(tmp)
        except Exception:
            pass
        return chrn


"""
if __name__ == "__main__":
    wt = "WT.out"
    z2 = WTOut(path=wt).get_z2_index()
    print(z2)
    chrn = WTOut(path=wt).get_chern_number()
    print(chrn)
"""
