import pandas as pd
import numpy as np



def get_parity(file, num_occ):
    vhkl = pd.read_csv("data/TIProj_candidates/YCr4Cu3O12/"+file+"_VHKL1.txt", header=None, delim_whitespace=True)
    # vhkl = pd.read_csv("data/TIProj_candidates/Si/1__VHKL.txt", header=None, delim_whitespace=True)
    pairs = []
    paired_index = []
    missing = []
    for index, row in vhkl.iterrows():
        vec = list(row)
        #if file == "6":
        #    vec[0] = vec[0]-0.5
        if vec[0] != 0 or vec[1] != 0 or vec[2] != 0:
            par = vhkl.index[(vhkl[0] == -vec[0]) & (vhkl[1] == -vec[1]) & (vhkl[2] == -vec[2])].tolist()
            dups = vhkl.index[(vhkl[0] == vec[0]) & (vhkl[1] == vec[1]) & (vhkl[2] == vec[2])].tolist()
            if len(par) > 1:
                print("more than one equal vector: " + str(par))
            elif len(dups) > 1:
                print("Duplicate vectors: " + str(dups) + " " + str(vec))
            elif len(par) == 0:
                missing.append(vec)
            elif index not in paired_index:
                pairs.append([index, par[0]])
                paired_index.append(index)
                paired_index.append(par[0])
    #print(len(paired_index))
    #print(len(vhkl))

    df_s1 = pd.read_csv("data/TIProj_candidates/YCr4Cu3O12/" + file + "_GWFS1.txt", header=None)
    df_s2 = pd.read_csv("data/TIProj_candidates/YCr4Cu3O12/" + file + "_GWFS2.txt", header=None)
    df_s1 = df_s1.iloc[:, :-1]
    df_s2 = df_s2.iloc[:, :-1]
    #print(missing)
    #print(df.shape)

    parity = []
    zero = complex("0+0j")
    for band in np.arange(0, num_occ, 2):
        coef_vals_s1 = list(df_s1.iloc[band])
        coef_vals_s2 = list(df_s2.iloc[band+1])
        vals = []
        standard = 0
        for pair in pairs:
            values = [complex(coef_vals_s1[pair[0]].replace(" ", "")), complex(coef_vals_s2[pair[1]].replace(" ", ""))]
            if not (values[0] == zero and values[1] == zero):
                if values[0] == zero or values[1] == zero:
                    print("err, one val zero")
                else:
                    #ratio = values[0] / values[1]
                    cutoff = 0.0
                    if np.sqrt(np.real(values[0])**2+np.imag(values[0])**2) > cutoff and np.sqrt(np.real(values[1])**2+np.imag(values[1])**2) > cutoff:
                        ratio = np.sqrt(np.real(values[0])**2+np.imag(values[0])**2) / np.sqrt(np.real(values[1])**2+np.imag(values[1])**2)
                        #print(ratio)
                        if np.real(values[0])*np.real(values[1]) < 0 and -1*np.imag(values[0])*np.imag(values[1]) < 0:
                            vals.append(round(-1*ratio,4))
                        elif np.real(values[0])*np.real(values[1]) > 0 and -1*np.imag(values[0])*np.imag(values[1]) > 0:
                            vals.append(round(ratio, 4))
                        else:
                            1#print("Problem with: " + str(values))
        #print(vals)
        parity.append(np.round(np.average(vals), 4))
    #print(len(parity))
    print(parity)
    print(np.prod(parity))
    return parity

if __name__ == "__main__":
    num_occ = 190
    par =  [get_parity("1",num_occ)]
    #par.append(get_parity("6",num_occ))
    #par.append(get_parity("47",num_occ))
    #par.append(get_parity("52",num_occ))
    #par.append(get_parity("277",num_occ))
    #par.append(get_parity("282", num_occ))

    for i, p1 in enumerate(par):
        #print(np.sum(np.array(p1)))
        for p2 in par[i+1:]:
            #print(np.array(p1) - np.array(p2))
            #print(np.sum(np.array(p1) - np.array(p2)))
            1+1

