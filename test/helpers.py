
def compare_files(f1,f2, skip_lines=0):
    with open(f1,"r") as newfile, open(f2, "r") as oldfile:
        for _ in range(skip_lines):
            newfile.readline()
            oldfile.readline()
        assert newfile.read() == oldfile.read()
