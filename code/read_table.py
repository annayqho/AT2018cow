""" Read the table from the radio data """


def zauderer():
    lines = open("table.txt", "r").readlines()
    dt = []
    nu = []
    f = []
    islim = []
    ef = []

    for line in lines:
        dt.append(float(line.split("&")[1]))
        nu.append(float(line.split("&")[3]))
        if '<' in line:
            islim.append(True)
            f.append(float(
                line.split("&")[4].split("$")[1][1:]))
        else:
            islim.append(False)
            f.append(float(
                line.split("&")[4].split("$")[1].split("pm")[0]))
            ef.append(float(
                line.split("&")[4].split("$")[1].split("pm")[1]))

    dt = np.array(dt)
    nu = np.array(nu)
    f = np.array(f)
    ef = np.array(ef)
    islim = np.array(islim)

    return nu, dt, f, ef, islim
