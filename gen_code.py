    # uint32_t depth = r;
    # uint32_t curr_size = size;
    # uint32_t curr_base = sum_prev_size;
    # uint32_t cbase2 = curr_base*2;
    # uint32_t *off_i = ij_buf + cbase2;
    # uint32_t *off_j = ij_buf + cbase2 + curr_size;

    # i_0 = off_i[curr_node]
    # j_0 = off_j[curr_node]

    # depth -= 1
    # curr_size = size_buf[depth];
    # curr_base -= curr_size;
    # uint32_t cbase2 = curr_base*2;
    # uint32_t *off_i = ij_buf + cbase2;
    # uint32_t *off_j = ij_buf + cbase2 + curr_size;

    # i_0_0 = off_i[i_0]
    # i_0_1 = off_j[i_0]
    # j_0_0 = off_i[j_0]
    # j_0_1 = off_j[j_0]

    # curr_size = size_buf[r-1];
    # r -= 1;
    # curr_base -= curr_size;

    # i_0_0_0 = *(ij_buf + curr_base*2 + i_0_0)
    # i_0_0_1 = *(ij_buf + curr_base*2 + curr_size + i_0_0)
    # i_0_1_0 = *(ij_buf + curr_base*2 + i_0_1)
    # i_0_1_1 = *(ij_buf + curr_base*2 + curr_size + i_0_1)
    # j_0_0_0 = *(ij_buf + curr_base*2 + j_0_0)
    # j_0_0_1 = *(ij_buf + curr_base*2 + curr_size + j_0_0)
    # j_0_1_0 = *(ij_buf + curr_base*2 + j_0_1)
    # j_0_1_1 = *(ij_buf + curr_base*2 + curr_size + j_0_1)

vals = ['i_0', 'j_0']
new_vals = []

print '''
uint32_t depth = r;
uint32_t curr_size = size;
uint32_t curr_base = sum_prev_size;
uint32_t cbase2 = curr_base*2;
uint32_t *off_i = ij_buf + cbase2;
uint32_t *off_j = ij_buf + cbase2 + curr_size;

i_0 = off_i[curr_node]
j_0 = off_j[curr_node]'''

rounds = range(1, 4)
for r in rounds:
    print '''
depth -= 1
curr_size = size_buf[depth];
curr_base -= curr_size;
uint32_t cbase2 = curr_base*2;
uint32_t *off_i = ij_buf + cbase2;
uint32_t *off_j = ij_buf + cbase2 + curr_size;
'''

    count = 0
    for ea in vals:
        if r == rounds[-1]:
            base0 = 'roots[%d] = '
            base1 = 'roots[%d] = '
            fmat0 = (count*2, ea)
            fmat1 = (count*2+1, ea)
        else:
            base0 = 'uint32_t %s_0 = '
            base1 = 'uint32_t %s_1 = '
            fmat0 = (ea, ea)
            fmat1 = (ea, ea)

        print (base0 + 'off_i[%s];') % fmat0
        print (base1 + 'off_j[%s];') % fmat1

        new_vals.append('%s_0' % ea)
        new_vals.append('%s_1' % ea)
        count += 1
    vals = new_vals
    new_vals = []
