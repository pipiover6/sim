def similar(opt1, opt2, turn):
    for i in range(3):
        if i != turn and opt1[i] != opt2[i]:
            return False
    return True

def subp(options, turn):
    ret = []
    for opt in options:
        similars = 0
        for opt2 in options:
            if similar(opt, opt2, turn):
                similars += 1
        assert(similars > 0)
        if similars == 1:
            ret.append(opt)
    return ret

def conway_paterson():
    sums = [6,7,8]
    options = []
    for i in range(9):
        for j in range(9):
            for k in range(9):
                if i + j + k in sums:
                    options.append((i,j,k))
    turn = 0
    while len(options) > 0:
        print(f"remaining at turn {turn}: {options}")
        knowable = subp(options, turn % 3)
        turn += 1
        for opt in knowable:
            print(f"{turn} {opt}")
            options.remove(opt)

conway_paterson()
