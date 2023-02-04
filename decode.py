#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 21:38:11 2022

@author: lukas
"""

import itertools
import matplotlib.pyplot as plt
import scipy.io
import numpy as np
#%%
def extract_subcode(frame) :
    return frame[24+3:24+3+14]

def list_to_int(l) :
    return sum(a*(1<<b) for a,b in zip(l, reversed(range(len(l)))))

def printbin(l) :
    print("".join("01"[x] for x in l))

def bcd_to_dec(x) :
    return (x&0xf) + 10*((x>>4)&0xf)

#%%

with open('data/frames.txt', 'r') as infile:
    frames = [[bool(int(x)) for x in line.strip()] for line in infile]
    
#%%

s0 = [x == "1" for x in '00100000000001']
s1 = [x == "1" for x in '00000000010010']

blocks = []

for i, frame in enumerate(frames) :
    subcode = extract_subcode(frame)
    if subcode == s0 :
        #print(i, "s0")
        if extract_subcode(frames[i+1]) == s1 :
            b = frames[i:i+98]
            if len(b) == 98 :
                blocks.append(b)
    elif subcode == s1 :
        #print(i, "s1")
        pass

#%%

efms = {}

with open("efm.txt") as infile :
    for line in infile:
        dec, b, efm = line.strip().split()
        dec = int(dec)
        b = int(b, 2)
        efm = int(efm, 2)
        #assert(dec == b)
        efms[efm] = b
        #print(line)
#%%

def calc_crc(bits) :
    poly = 0x1021
    crc = 0
    for bit in bits :
        if (crc>>15)&1 != bit :
            crc = ((crc<<1)&0xffff) ^ poly
        else :
            crc = ((crc<<1)&0xffff)
    return crc
    
    

for block in blocks :
    subcodes = [extract_subcode(frame) for frame in block]
    assert(subcodes[0] == s0)
    assert(subcodes[1] == s1)
    subcodes_payload = [efms.get(list_to_int(x), 0) for x in subcodes[2:]]
    assert(len(subcodes_payload) == 96)
    q = [bool(x & 64) for x in subcodes_payload]
    q_for_crc = [ not x if i>= 96-16 else x for i,x in enumerate(q)]
    #printbin(q_for_crc)
    crc_syndrome = calc_crc(q_for_crc)
    adr = list_to_int(q[4:4+4])
    #print(adr)
    if crc_syndrome == 0 :
        if adr == 1:
            data_bytes = [list_to_int(q[x:x+8]) for x in range(8, 80, 8)]       
            tno, idx, rmin, rsec, rframe, zero, amin, asec, aframe = (bcd_to_dec(x) for x in data_bytes)
            
            print(f"Track {tno}.{idx} R={rmin:02d}:{rsec:02d}:{rframe:02d} A={amin:02d}:{asec:02d}:{aframe:02d}")
        else :
            print("A", adr)
    else :
        print("CRC error")
    #printbin(q)
    
#%%


def symbols_from_frame(frame) :
    symbol_idxs = (24+3+14+3+(i*(14+3)) for i in range(32))
    return [efms[list_to_int(frame[i:i+14])] for i in symbol_idxs]

class Delay:
    def __init__(self, n_delay, fill = None) :
        self.register = [fill]*n_delay
        
    def step(self, v) :
        if len(self.register) == 0 :
            return v
        r = self.register[-1]
        self.register = [v] + self.register[:-1]
        return r
        

deinterleave_tab = (
0,
1,
6,
7,
16,
17,
22,
23,
2,
3,
8,
9,
18,
19,
24,
25,
4,
5,
10,
11,
20,
21,
26,
27,
)

def has_last_delay(i) :
    return i in (4,5,6,7, 12,13,14,15, 20,21,22,23)

def extract_audio(frames) :
    first_delays = [Delay(0 if i%2 == 0 else 1, 0) for i in range(32)]
    second_delays = [Delay(i*4, 0) for i in reversed(range(28))]
    third_delays = [Delay(2 if has_last_delay(i) else 0, 0) for i in range(24)]
    for frameidx, frame in enumerate(frames) :
        try :
            symbols_in = symbols_from_frame(frame)

            symbols_delayed1 = [d.step(x) for d,x in zip(first_delays, symbols_in)]
            #symbols_delayed1_inverted = [x != (i in (12,13,14,15, 28,29,30,31)) for i,x in enumerate(symbols_delayed1)]
            
            # skip c1 decoder
            
            symbols_delayed2 = [d.step(x) for d,x in zip(second_delays, symbols_delayed1)]
            
            #skip c2 decoder
            
            symbols_deinterleaved = [symbols_delayed2[i] for i in deinterleave_tab]
            
            symbols_out = [d.step(x) for d,x in zip(third_delays, symbols_deinterleaved)]
            
            yield from symbols_out
            
        except KeyError as e :
            print(i,e)


def combine_samples(samples_raw) :
    while chunk := list(itertools.islice(samples_raw, 2)):
        if len(chunk) != 2 :
            continue
        i = (chunk[0] << 8) + chunk[1]
        if i & (1<<15) :
            i -= (1<<16)
        yield i


#%%
samples = np.array(list(combine_samples(extract_audio(frames))), dtype=np.int16)
scipy.io.wavfile.write("data/cd.wav", 44100, samples.reshape((len(samples)//2, 2)))
