def name_extracter(file, typo):
    print(file.readline().strip())
    if typo == 'int':
        var = int(file.readline().strip())
    elif typo == 'flo':
        var = float(file.readline().strip())
    elif typo == 'str':
        var = str(file.readline().strip())
    print(var)
    return var

def input_reader():
    with open('./input.txt', 'r') as f:
        print('----------INPUT DATA----------')
        flamenum = name_extracter(f, 'int')
        zst = name_extracter(f, 'flo')
        press = name_extracter(f, 'flo')
        tempf = name_extracter(f, 'flo')
        tempo = name_extracter(f, 'flo')
        mdot = name_extracter(f, 'flo')
        compf = name_extracter(f, 'str')
        compo = name_extracter(f, 'str')
        mech = name_extracter(f, 'str')
        outname = name_extracter(f, 'str')
        moleout = name_extracter(f, 'str')
        width = name_extracter(f, 'flo')
        npoint = name_extracter(f, 'int')
        nzeta = name_extracter(f, 'int')
        print(f.readline().strip())
        cdef = name_extracter(f, 'str').split(',')
        print('----------END OF DATA----------')
    print('Is this OK?[y/n]')
    usrinput = input()
    if usrinput == 'n':
        print('Cancelling...')
    else:
        input_package = (flamenum, zst, press, tempf, tempo, mdot, compf, compo, mech, outname, moleout, width, npoint, nzeta, cdef)
        return input_package
