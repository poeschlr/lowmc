SBox = [0x00, 0x01, 0x03, 0x06, 0x07, 0x04, 0x05, 0x02]
sbox_size = 3

ddt = [[0 for j in range(2**sbox_size)] for i in range(2**sbox_size)]

def print_list(list, name):
    print("{:s} = [".format(name))
    for i in range(len(list)):
        print("[{:s}]{:s}".format(', '.join(map(str,list[i])),',' if i<len(list)-1 else ']'))


#print_ddt()

for i1 in range(2**sbox_size):
    for id in range(2**sbox_size):
        od = SBox[i1]^SBox[i1^id]
        ddt[id][od] += 1



possible_out_d = []
possible_in_d = []


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

print_list(ddt, 'ddt')
print_list(possible_out_d, 'possible_out_d')
print_list(possible_in_d, 'possible_in_d')
