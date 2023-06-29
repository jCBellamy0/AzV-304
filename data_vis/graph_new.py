import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.transforms as tran
import sys
import collections
class FullDataset:
    def __init__(self, filename, eff_temp=-1.0, log_g=-1.0):
                self.eff_temp = eff_temp
                self.log_g = log_g
                self.storage, self.error_storage = self.read_data(filename)
    def read_data(self, filename):#stored in dict, first num indicates the micro/abund combo, second indicates ionization stage, third is the line location
        f = open(filename, "r")
        f_lines = f.readlines()
        f.close()
        storage = dict()
        error_storage=dict()
        for i in range(27):
                storage[i]=dict()
                error_storage[i]=dict()
        # 12-Line location, 4-EW, 13-Transition 
        read_cols = [13, 12, 4]
        count=0
        for l in f_lines[30:]:
            line = l.split()
            # if labeled species
            if len(line) > 13 and '(' in line[13]:
                    k = line[13][:line[13].index('(')]
                    if not k in storage[count%27]:
                            storage[count%27][k] = dict()
                            error_storage[count%27][k] = dict()
                    storage[count%27][k][float(line[12])]=float(line[4])
                    error_storage[count%27][k][float(line[12])]=float(line[5])
            count+=1
        return storage, error_storage
    def get_micro_abund_storage(self, micro, abund):#returns sub-dict of micro/abund combo
        abund=str(abund)
        abund_convert={"":0,"C/2":1,"C*2":2,"N/2":3,"N*2":4,"O/2":5,"O*2":6,"SI/2":7,"SI*2":8}
        if micro == 1:
            if abund.isdigit():
                return self.storage[0+int(abund)]
            else:
                return self.storage[0+abund_convert[abund]]
        elif micro == 3:
            if abund.isdigit():
                return self.storage[9+int(abund)]
            else:
                return self.storage[9+abund_convert[abund]]
        elif micro == 5:
            if abund.isdigit():
                return self.storage[18+int(abund)]
            else:
                return self.storage[18+abund_convert[abund]]
        else:
                print("Error: incorrect microturbulence value entered")
                sys.exit(0)
    def get_micro_abund_error_storage(self, micro, abund):#returns sub-dict of micro/abund combo
        abund=str(abund)
        abund_convert={"":0,"C/2":1,"C*2":2,"N/2":3,"N*2":4,"O/2":5,"O*2":6,"SI/2":7,"SI*2":8}
        if micro == 1:
            if abund.isdigit():
                return self.error_storage[0+int(abund)]
            else:
                return self.error_storage[0+abund_convert[abund]]
        elif micro == 3:
            if abund.isdigit():
                return self.error_storage[9+int(abund)]
            else:
                return self.error_storage[9+abund_convert[abund]]
        elif micro == 5:
            if abund.isdigit():
                return self.error_storage[18+int(abund)]
            else:
                return self.error_storage[18+abund_convert[abund]]
        else:
                print("Error: incorrect microturbulence value entered")
                sys.exit(0)


class ObsDataset:
    def __init__(self, filename):
            self.storage, self.error_storage = self.read_data(filename)
    def read_data(self, filename):
            f = open(filename, "r")
            f_lines = f.readlines()
            f.close()
            storage = dict()
            error_storage = dict()
            # 12-Line location, 4-EW, 13-Transition
            read_cols = [13, 12, 4]
            for l in f_lines[4:]:
                    line = l.split()
                    # if labeled species
                    if len(line) > 13 and '(' in line[13]:
                            k = line[13][:line[13].index('(')]
                            if not k in storage:
                                    storage[k] = dict()
                                    error_storage[k]=dict()
                            storage[k][float(line[12])]=float(line[4])
                            error_storage[k][float(line[12])]=float(line[5])
            return storage, error_storage

def get_ionization_counts(obs_dataset, other_dataset):
    count_dict={}    
    for species in obs_dataset:
            if species in other_dataset:
                if species not in count_dict:
                    count_dict[species]=0
                for line in obs_dataset[species]:
                    if line in other_dataset[species]:
                        count_dict[species]=count_dict[species]+1
    return count_dict


def get_chi_val(obs_dataset, other_dataset):
    sum = 0.0
    count_dict = get_ionization_counts(obs_dataset.storage,other_dataset)
    for species in obs_dataset.storage:
            if species in other_dataset and species != "Ca2":
                for line in obs_dataset.storage[species]:
                    if line in other_dataset[species]:
                        ew1 = obs_dataset.storage[species][line]
                        ew2 = other_dataset[species][line]
                        error_val = obs_dataset.error_storage[species][line]
                        sum += float(float((ew1-ew2))/error_val)**2
                        #sum += float(float(float((ew1-ew2)**2)/error_val)/count_species
    return round(sum, 2)


def get_chi_val_interp(obs_dataset, other_dataset1, other_dataset2):
    sum = 0.0
    count_dict = get_ionization_counts(obs_dataset.storage,other_dataset1)
    for species in obs_dataset.storage:
        if species in other_dataset1 and species in other_dataset2 and species != "Ca2":
            for line in obs_dataset.storage[species]:
                if line in other_dataset1[species] and line in other_dataset2[species]:
                    ew1 = obs_dataset.storage[species][line]
                    ew2 = other_dataset1[species][line]
                    ew3 = other_dataset2[species][line]
                    ex = float((ew2+ew3)/2)
                    error_val = obs_dataset.error_storage[species][line]
                    sum += float(float((ew1-ex))/error_val)**2

    return round(sum, 2)

def plot_chi_wavelength(obs_dataset, other_dataset,filter_list=None):#obs_storage, other_dataset.storage list
        pos_x_vals = []
        pos_y_vals = []
        pos_species_line = []
        neg_x_vals = []
        neg_y_vals = []
        neg_species_line = []

        count_dict = get_ionization_counts(obs_dataset,other_dataset)
        print(count_dict)
        for species in obs_dataset:
            if filter_list == None and species in other_dataset:
                for line in obs_dataset[species]:
                    if line in other_dataset[species]:
                        ew1 = obs_dataset[species][line]
                        ew2 = other_dataset[species][line]
                        if (ew1-ew2)>=0.0:
                            neg_x_vals.append(line)
                            neg_y_vals.append(float(float(float((ew1-ew2)**2)/ew2)/count_dict[species]))
                            neg_species_line.append(species+" "+str(line))
                        else:
                            pos_x_vals.append(line)
                            pos_y_vals.append(float(float(float((ew1-ew2)**2)/ew2)/count_dict[species]))
                            pos_species_line.append(species+" "+str(line))
            elif species in filter_list and species in other_dataset:
                for line in obs_dataset[species]:
                    if line in other_dataset[species]:
                        ew1 = obs_dataset[species][line]
                        ew2 = other_dataset[species][line]
                        if (ew1-ew2)>=0.0:
                            neg_x_vals.append(line)
                            neg_y_vals.append(float(float(float((ew1-ew2)**2)/ew2)/count_dict[species]))
                            neg_species_line.append(species+" "+str(line))
                        else:
                            pos_x_vals.append(line)
                            pos_y_vals.append(float(float(float((ew1-ew2)**2)/ew2)/count_dict[species]))
                            pos_species_line.append(species+" "+str(line))
            else:
                continue


        plt.scatter(pos_x_vals,pos_y_vals,c='blue',label="Overfit")
        for i, txt in enumerate(pos_species_line):
            plt.annotate(txt,(pos_x_vals[i],pos_y_vals[i]),xytext=(0.0,65.0),textcoords='offset points',rotation=90)
        plt.scatter(neg_x_vals,neg_y_vals,c='red',label="Underfit")
        for i, txt in enumerate(neg_species_line):
            plt.annotate(txt,(neg_x_vals[i],neg_y_vals[i]),xytext=(0.0,65.0),textcoords='offset points',rotation=90)
        plt.axhline(y = 0.0, color = 'black', linestyle = 'dashed')
        plt.ylim([-1,5])
        plt.xlabel("Wavelength (Angstroms)")
        plt.ylabel("Chi-square contribution")
        plt.legend(loc ="upper right")
        plt.show()

def plot_error_ew(obs_dataset):#obs_storage, other_dataset.storage list 
        x_vals = []
        y_vals = []

        #print(count_dict)
        ews = obs_dataset.storage
        errors = obs_dataset.error_storage
        for species in ews:
                for line in ews[species]:
                        x_vals.append(ews[species][line])
                        y_vals.append(errors[species][line])
        plt.scatter(x_vals,y_vals,c='blue')
        plt.xlabel("EW (mA)")
        plt.ylabel("Error (%)")
        plt.xscale("log")
        plt.show()


obs_dataset = ObsDataset("EW_DATA_REDUX2")

mod_dataset_1 = FullDataset("EW_DATA_FULL_3p52p75")
mod_dataset_2 = FullDataset("EW_DATA_FULL_3p62p75")
mod_dataset_3 = FullDataset("EW_DATA_FULL_3p52p7")
mod_dataset_4 = FullDataset("EW_DATA_FULL_3p62p7")
mod_dataset_5 = FullDataset("EW_DATA_FULL_3p72p7")
mod_dataset_6 = FullDataset("EW_DATA_FULL_3p52p6")
mod_dataset_7 = FullDataset("EW_DATA_FULL_3p62p6")
mod_dataset_8 = FullDataset("EW_DATA_FULL_3p72p6")

def get_chi_val_abund(obs_dataset, mod_dataset, mt, abund):
    sum = 0.0
    other_dataset = mod_dataset.get_micro_abund_storage(mt,abund)
    count_dict = get_ionization_counts(obs_dataset.storage,other_dataset)
    for species in obs_dataset.storage:
            if species in other_dataset and (species == "C2" or species == "CIII" or species == "CIV"):
                carbon_dataset = mod_dataset.get_micro_abund_storage(mt,"C*2")
                for line in obs_dataset.storage[species]:
                    if line in carbon_dataset[species]:
                        ew1 = obs_dataset.storage[species][line]
                        ew2 = carbon_dataset[species][line]
                        error_val = obs_dataset.error_storage[species][line]
                        sum += float(float((ew1-ew2))/error_val)**2
            elif species in other_dataset and species != "Ca2":
                for line in obs_dataset.storage[species]:
                    if line in other_dataset[species]:
                        ew1 = obs_dataset.storage[species][line]
                        ew2 = other_dataset[species][line]
                        error_val = obs_dataset.error_storage[species][line]
                        sum += float(float((ew1-ew2))/error_val)**2            

    return round(sum, 2)


def get_chi_val_abund_interp(obs_dataset, mod_dataset, mt, abund):
    sum = 0.0
    other_dataset_1 = mod_dataset.get_micro_abund_storage(mt-1,abund)
    other_dataset_2 = mod_dataset.get_micro_abund_storage(mt+1,abund)
    count_dict = get_ionization_counts(obs_dataset.storage,other_dataset_1)
    for species in obs_dataset.storage:
           if species in other_dataset_1 and species in other_dataset_2 and (species == "C2" or species == "CIII" or species == "CIV"):
                carbon_dataset_1 = mod_dataset.get_micro_abund_storage(mt-1,"C*2")
                carbon_dataset_2 = mod_dataset.get_micro_abund_storage(mt+1,"C*2")
                for line in obs_dataset.storage[species]:
                    if line in carbon_dataset_1[species] and line in carbon_dataset_2[species]:
                        ew1 = obs_dataset.storage[species][line]
                        ew2 = carbon_dataset_1[species][line]
                        ew3 = carbon_dataset_2[species][line]
                        ex = float((ew2+ew3)/2)
                        error_val = obs_dataset.error_storage[species][line]
                        sum += float(float((ew1-ex))/error_val)**2
           if species in other_dataset_1 and species in other_dataset_2 and species != "Ca2":
                for line in obs_dataset.storage[species]:
                    if line in other_dataset_1[species] and line in other_dataset_2[species]:
                        ew1 = obs_dataset.storage[species][line]
                        ew2 = other_dataset_1[species][line]
                        ew3 = other_dataset_2[species][line]
                        ex = float((ew2+ew3)/2)
                        error_val = obs_dataset.error_storage[species][line]
                        sum += float(float((ew1-ex))/error_val)**2
    return round(sum, 2)

def plot_best_fit_abund_chi_value(obs_dataset, mod_dataset):
    txt_table=[[],[],[],[],[]]#5 mt rows, 9 abund cols
    for i in range(1,6):
        if i in [1,3,5]:
            for j in range(9):
                txt_table[i-1].append(get_chi_val_abund(obs_dataset,mod_dataset,i,j))
        else:
            for j in range(9):
                txt_table[i-1].append(get_chi_val_abund_interp(obs_dataset,mod_dataset,i,j))

    colors=[[],[],[],[],[]]
    lowest_val = 10000000000000
    lowest_indexes = (0,0)
    for mt in range(5):
        for abund in range(9):
            if(txt_table[mt][abund]<lowest_val):
                lowest_val=txt_table[mt][abund]
                lowest_indexes=(mt,abund)
            colors[mt].append("white")
    colors[lowest_indexes[0]][lowest_indexes[1]]='green'
    return (txt_table, colors)


def plot_chi_square(obs_dataset,mod_dataset):
    txt_table=[[],[],[],[],[]]#5 mt rows, 9 abund cols
    for i in range(1,6):
        if i in [1,3,5]:
            for j in range(9):
                txt_table[i-1].append(get_chi_val(obs_dataset,mod_dataset.get_micro_abund_storage(i,j)))
        else:
            for j in range(9):
                txt_table[i-1].append(get_chi_val_interp(obs_dataset,mod_dataset.get_micro_abund_storage(i-1,j),mod_dataset.get_micro_abund_storage(i+1,j)))

    colors=[[],[],[],[],[]]
    lowest_val = 10000000000000
    lowest_indexes = (0,0)
    for mt in range(5):
        for abund in range(9):
            if(txt_table[mt][abund]<lowest_val):
                lowest_val=txt_table[mt][abund]
                lowest_indexes=(mt,abund)
            colors[mt].append("white")
    colors[lowest_indexes[0]][lowest_indexes[1]]='green'
    return (txt_table, colors)

def plot_chi_table(obs_dataset, table_func, mod_dataset_list):
        fig = plt.figure(1)
        ax_main = fig.add_subplot(111) 
        count=1
        for mod_dataset in mod_dataset_list:
            if count==3:
                count=4
            ax = fig.add_subplot(3,3,count)
            ax.axis("off")
            ax.axis("tight")
            ax_data = table_func(obs_dataset, mod_dataset)
            tbl = ax.table(cellText=ax_data[0],cellColours=ax_data[1],rowLabels=["MT1", "MT2", "MT3", "MT4", "MT5"], colLabels=["N/A","C/2","C*2","N/2","N*2","O/2","O*2","SI/2","SI*2"],loc='center')
            font_size=0.5
            tbl.auto_set_font_size(False)
            #tbl.set_fontsize(font_size)
            #tbl.scale(1.7,2.0)
            count+=1
        ax_main.spines['top'].set_color('none')
        ax_main.spines['bottom'].set_color('none')
        ax_main.spines['left'].set_color('none')
        ax_main.spines['right'].set_color('none')
        ax_main.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
        ax_main.set_xlabel('Log G: 3.5 - 3.7')
        ax_main.set_ylabel('T_Eff: 2.6 - 2.75')
        #fig.patch.set_visible(False)
        plt.show()

def check_ew_comparison(mod_full_dataset):
     helper_dict = {1:["C2","CIII","CIV"],3:["NIII"],5:["OIII","O2","OIV"],7:["SkIII","SkIV"]}
     for mt_val in range(1,7,2):
        for abund_val in range(1,8,2):
            other_dataset_norm = mod_full_dataset.get_micro_abund_error_storage(mt_val,0)
            other_dataset_half = mod_full_dataset.get_micro_abund_error_storage(mt_val,abund_val)
            other_dataset_double = mod_full_dataset.get_micro_abund_error_storage(mt_val,abund_val+1)
            for species in other_dataset_norm:
                if species in other_dataset_half and species in other_dataset_double:
                    for line in other_dataset_norm[species]:
                        if line in other_dataset_half[species] and line in other_dataset_double[species]:
                            ew_norm = other_dataset_norm[species][line]
                            ew_half = other_dataset_half[species][line]
                            ew_double = other_dataset_double[species][line]
                            if (species in helper_dict[abund_val]) and (ew_half>ew_norm or ew_norm>ew_double):
                                print("MT: "+str(mt_val)+" Species: "+str(species)+" Line: "+str(line))
                                print("EW/2: "+str(ew_half)+" EW: "+str(ew_norm)+" EW*2: "+str(ew_double))


#plot_chi_wavelength(obs_dataset.storage,mod_dataset_4.get_micro_abund_storage(1,"C*2"))

#plot_error_ew(obs_dataset)

#print(get_ionization_counts(obs_dataset.storage,mod_dataset_1.storage[0]))

plot_chi_table(obs_dataset, plot_chi_square, [mod_dataset_1, mod_dataset_2,mod_dataset_3,mod_dataset_4,mod_dataset_5,mod_dataset_6,mod_dataset_7,mod_dataset_8])

#check_ew_comparison(mod_dataset_4)
