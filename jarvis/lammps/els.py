import glob, json

mem = []
for file in glob.glob("*.alloy"):
    ff = str(file)
    element_ff = []
    f = open(ff, "r")
    list_el = []
    lines = f.readlines()
    content = (lines[3]).split(" ")
    # content=(lines[3]).split("' '|\n|\r\n")
    for val in content:

        if val != "" and val != "\n" and val != "\r\n":
            list_el.append(val)
    for i in range(0, len(list_el)):
        if i != 0:
            # element_ff.append(list_el[i])
            element_ff.append(list_el[i].split("\n")[0])
    print(file, element_ff)
    info = {}
    info["elements"] = element_ff
    info["file"] = file
    mem.append(info)
print(len(mem))
f = open("alloy_pots.json", "w")
f.write(json.dumps(mem, indent=4))
f.close()
