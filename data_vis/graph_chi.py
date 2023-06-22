import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.transforms as tran
import sys
# read in obs,mt1_on2, mt1_tim2,mt3,mt5


# returns mt vel and abundance as int, str based on filename structure "EW_DATA_MT[VALUE]_[SPECIES]_[ADJUSTMENT]"
def parse_filename(filename):
    if len(filename) < 11:
        return -1, ""
    elif len(filename) < 12:
        return int(filename[10]), "N/A"
    return int(filename[10]), filename[12:]


def get_chi_val(obs_dataset, other_dataset):
    sum = 0.0
    for species in obs_dataset.storage:
            if species in other_dataset.storage:
                for line in obs_dataset.storage[species]:
                    if line in other_dataset.storage[species]:
                        ew1 = obs_dataset.storage[species][line]
                        ew2 = other_dataset.storage[species][line]
                        sum += float(float((ew1-ew2)**2)/ew2)
    return round(sum, 3)


def get_chi_val_interp(obs_dataset, other_dataset1, other_dataset2):
    sum = 0.0
    for species in obs_dataset.storage:
        if species in other_dataset1.storage and species in other_dataset2.storage:
            for line in obs_dataset.storage[species]:
                if line in other_dataset1.storage[species] and line in other_dataset2.storage[species]:
                    ew1 = obs_dataset.storage[species][line]
                    ew2 = other_dataset1.storage[species][line]
                    ew3 = other_dataset2.storage[species][line]
                    ex = float((ew2+ew3)/2)
                    sum += float(float((ew1-ex)**2)/ex)
    return round(sum, 3)
# reads LINE_ID for the optical depth, returns species-sorted dict of dicts of wavelength


def read_depth_data(filename):
    f = open(filename, "r")
    f_lines = f.readlines()
    f.close()
    storage = dict()
    read_cols = [6, 1, 2]
    for l in f_lines:
        line = l.split()
        if '(' in line[6]:
            k = line[6][:line[6].index('(')]
            if not k in storage:
                storage[k] = dict()
            storage[k][float(line[1])] = float(line[2])
    return storage



class Dataset:
    def __init__(self, filename, is_obs, eff_temp=-1.0, log_g=-1.0):
                if (is_obs):
                        self.eff_temp = -1.0
                        self.log_g = -1.0
                        self.mt_vel = -1
                        self.abund = ""
                        self.storage = self.rep_ids(filename)
                else:
                        self.eff_temp = eff_temp
                        self.log_g = log_g
                        self.mt_vel, self.abund = parse_filename(filename)
                        self.storage = self.read_data(filename)

    def read_data(self, filename):
                f = open(filename, "r")
                f_lines = f.readlines()
                f.close()
                storage = dict()
                # 12-Line location, 4-EW, 13-Transition
                read_cols = [13, 12, 4]
                for l in f_lines[4:]:
                        line = l.split()
                        # if labeled species
                        if len(line) > 13 and '(' in line[13]:
                                k = line[13][:line[13].index('(')]
                                if not k in storage:
                                        storage[k] = dict()
                                storage[k][float(line[12])]=float(line[4])
                return storage
            
    def read_line_id(self,filename):
        f = open(filename, "r")
        f_lines = f.readlines()
        f.close()
        storage = dict()
        read_cols = [1,10]
        for l in f_lines[1:]:
            line = l.split()
            storage[float(line[1])] = line[10]
        return storage

    def rep_ids(self,filename):#replace
        count=0
        line_store=self.read_line_id("LINE_ID")
        f = open(filename, "r")
        f_lines = f.readlines()
        f.close()

        storage = dict()
        for l in f_lines[4:]:
            txt_line = l.split()
            loc=float(txt_line[2])
            lowest_diff=1.0
            lowest_key=1.0
            for k in line_store.keys():
                if abs(loc-k)<0.1:
                    if lowest_diff>abs(loc-k):
                        lowest_diff=abs(loc-k)
                        lowest_key=k
            if lowest_key!=1.0:
                i=line_store[lowest_key][:line_store[lowest_key].index('(')]
                if not i in storage:
                    storage[i] = dict()
                storage[i][lowest_key]=float(txt_line[4])
                count+=1
        print(count)
        return storage

EFF_TEMP=2.7
LOG_G=3.6
add_ins=['_SI', '_O', '_C', '_N']
files=['EW_DATA_REDUX', 'EW_DATA_MT1_ON2', 'EW_DATA_MT1_TIM2',
    'EW_DATA_MT3_ON2', 'EW_DATA_MT3_TIM2', 'EW_DATA_MT5_ON2', 'EW_DATA_MT5_TIM2']
add_files=['EW_DATA_MT1', 'EW_DATA_MT3', 'EW_DATA_MT5']
data_storage=dict()  # dict sorted by abundance, then mt
obs_data=Dataset(files[0], True)
for f in add_files:  # reads in spectrums that don't have abundance changes
    d=Dataset(f, False, EFF_TEMP, LOG_G)
    if not d.abund in data_storage:
        data_storage[d.abund]=dict()
    data_storage[d.abund][d.mt_vel]=d
for f in files[1:]:  # reads in all spectra into storage
    for a in add_ins:
        f_str=f[:11]+a+f[11:]
        #print(f_str)
        d=Dataset(f_str, False, EFF_TEMP, LOG_G)
        if not d.abund in data_storage:
            data_storage[d.abund]=dict()
        data_storage[d.abund][d.mt_vel]=d

fig, ax=plt.subplots()

# hide axes
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')
#    mt 1-5km/s
# population
# no change
# /2,*2 for each
txt_table=[[], [], [], [], []]
colors=[[], [], [], [], []]
cols=[]
lowest_chi_val=1000000.0
lowest_location=(0,0)
print(obs_data.storage)
for abund in data_storage:
    cols.append(abund)
    for mt in data_storage[abund]:
        if mt > 1:
            chi_val=get_chi_val_interp(obs_data,data_storage[abund][mt-2], data_storage[abund][mt])
            txt_table[mt-2].append(chi_val)
            colors[mt-2].append('white')
            if chi_val<lowest_chi_val:
                lowest_chi_val=chi_val
                lowest_location=(mt-2,len(txt_table[mt-2])-1)
        chi_val=get_chi_val(obs_data, data_storage[abund][mt])
        txt_table[mt-1].append(chi_val)
        colors[mt-1].append('white')
        if chi_val<lowest_chi_val:
            lowest_chi_val=chi_val
            lowest_location=(mt-1,len(txt_table[mt-1])-1)

colors[lowest_location[0]][lowest_location[1]]='green'
table=plt.table(cellText=txt_table, cellColours=colors, rowLabels=["MT1", "MT2", "MT3", "MT4", "MT5"], colLabels=cols,loc='center')
plt.title("Effective Temp="+str(EFF_TEMP)+" Log G="+str(LOG_G))
plt.show()
