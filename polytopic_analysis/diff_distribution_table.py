SBox = [0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02]
invSBox = [0x00, 0x01, 0x07, 0x02, 0x05, 0x06, 0x03, 0x04]
sbox_size = 3

ddt = [[0 for j in range(2**sbox_size)] for i in range(2**sbox_size)]

iddt = [[0 for j in range(2**sbox_size)] for i in range(2**sbox_size)]

def print_list(list, name):
    print("{:s} = [".format(name))
    for i in range(len(list)):
        print("[{:s}]{:s}".format(', '.join(map(str,list[i])),',' if i<len(list)-1 else ']'))


#print_ddt()

for i1 in range(2**sbox_size):
    for id in range(2**sbox_size):
        od = SBox[i1]^SBox[i1^id]
        ddt[id][od] += 1

for i1 in range(2**sbox_size):
    for id in range(2**sbox_size):
        od = invSBox[i1]^invSBox[i1^id]
        iddt[id][od] += 1

possible_out_d = []
possible_in_d = []
possible_out2_d = []


for id in range(2**sbox_size):
    possible_out_d.append([])
    for od in range(2**sbox_size):
        if ddt[id][od] > 0:
            possible_out_d[id].append(od)

for od in range(2**sbox_size):
    possible_in_d.append([])
    for id in range(2**sbox_size):
        if ddt[id][od] > 0:
            possible_in_d[od].append(id)

for id in range(2**sbox_size):
    possible_out2_d.append([])
    for od in range(2**sbox_size):
        if iddt[id][od] > 0:
            possible_out2_d[id].append(od)

print_list(ddt, 'ddt')
print_list(possible_in_d, 'possible_in_d')
print_list(possible_out_d, 'possible_out_d')
print_list(possible_out2_d, 'possible_out2_d')

out_d_f = []
out_d_r = []
for id in range(2**sbox_size):
    out_d_f.append([])
    out_d_r.append([])
    for anchor in range(2**sbox_size):
        out_d_f[id].append(SBox[anchor^id]^SBox[anchor])
        out_d_r[id].append(invSBox[anchor^id]^invSBox[anchor])


print_list(out_d_f, 'ddiff_prob_table_forward')
print_list(out_d_r, 'ddiff_prob_table_backward')
