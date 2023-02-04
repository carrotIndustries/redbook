import itertools
import matplotlib.pyplot as plt
import numpy as np

def read_csv(filename) :
    for line in open(filename, "r") :
        if line.startswith("#") :
            continue
        yield float(line.strip().split(",")[1])

def interpolate(samples, factor) :
    last = next(samples)
    yield last
    for sample in samples :
        for i in range(factor) :
            t = (i+1)/factor
            yield last*(1-t) + sample*(t)
        last = sample

def slice_bits(samples) :
    threshold = 0
    hyst = .01
    for sample in samples :
        if sample > threshold:
            yield True
            threshold = -hyst
        else :
            yield False
            threshold = hyst


#%%
samples = list(itertools.islice(read_csv("data/scope.csv"), 200))
plt.plot(samples)
plt.xlabel("sample")
plt.ylabel("Voltage [V]")
plt.xlim(0, 200)
plt.grid()

#%%

samples_interpolated =  interpolate( read_csv("data/scope.csv"), 20)

all_bits = slice_bits(samples_interpolated)

acc = 0
acc_size = 1000
ftw = 0
ftw0 = 42.7
lastbit = False

bits = []
accs = []
sampled_bits = []
phase_deltas_filtered = []
phase_delta = 0
phase_delta_filtered = 0
integ = 0
decoded_bits = []
last_acc = 0
ftws = []
integs = []
integ_max = 927681

alpha = .005
Ki = .0000004
Kp = .001


debug = False

for bit in all_bits:
    if debug :
        bits.append(bit)
    if bit != lastbit : # input transition
        phase_delta = (acc_size/2 - acc)
    if acc < last_acc : # phase accumulator has wrapped around
        sampled_bits.append(bit)
    last_acc = acc
    phase_delta_filtered = phase_delta*alpha + phase_delta_filtered*(1-alpha)
    integ += phase_delta_filtered
    if integ > integ_max :
        integ = integ_max
    elif integ < -integ_max :
        integ = -integ_max
    integs.append(integ)
    
    ftw = ftw0 + phase_delta_filtered*Kp + integ * Ki
    if debug :
        phase_deltas_filtered.append(phase_delta_filtered/acc_size)
    lastbit = bit
    acc = (acc+ftw)%acc_size
    if debug and len(bits) > 200000 :
        break
    if debug :
        accs.append(acc/acc_size)
        ftws.append(ftw)
    
#%%
plt.rcParams["figure.figsize"] = (15,4)
if debug :
    plt.plot(np.array(phase_deltas_filtered[::10000])*360)
    plt.xlabel("sample")
    plt.ylabel("filtered phase error [degree]")
    plt.grid()
print(ftw)
#%%

if debug :
    plt.xlim(len(bits)-400, len(bits))
    plt.plot(accs, label="VCO phase (normalized)")
    plt.plot([0, len(bits)], [.5, .5], label="180°")
    plt.plot(bits, label="bits")
    plt.grid()
    plt.xlabel("sample")
    plt.legend()

#%%

nrz_bits = [a != b for a,b in zip(sampled_bits[1:], sampled_bits)]

syncpat = [x == "1" for x in "1000-0000-0001-0000-0000-0010".replace("-", "")]
last_i = None
frames = []
frame_len = 588
for i in range(len(nrz_bits)-len(syncpat)) :
    if nrz_bits[i:i+len(syncpat)] == syncpat :
        f = nrz_bits[i:i+frame_len]
        if len(f) == frame_len :
            frames.append(f)    
        if last_i is not None :
           if i-last_i != frame_len :
                print("short fraeme", i-last_i)
        last_i = i

#%%

with open("data/frames.txt", "w") as ofile:
    ofile.write("\n".join("".join("01"[x] for x  in frame) for frame in frames))
