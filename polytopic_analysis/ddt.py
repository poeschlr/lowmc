ddt = [
[8, 0, 0, 0, 0, 0, 0, 0],
[0, 2, 0, 2, 0, 2, 0, 2],
[0, 0, 2, 2, 0, 0, 2, 2],
[0, 2, 2, 0, 0, 2, 2, 0],
[0, 0, 0, 0, 2, 2, 2, 2],
[0, 2, 0, 2, 2, 0, 2, 0],
[0, 0, 2, 2, 2, 2, 0, 0],
[0, 2, 2, 0, 2, 0, 0, 2]]
possible_out_d = [
[0],
[1, 3, 5, 7],
[2, 3, 6, 7],
[1, 2, 5, 6],
[4, 5, 6, 7],
[1, 3, 4, 6],
[2, 3, 4, 5],
[1, 2, 4, 7]]
possible_in_d = [
[0],
[1, 3, 5, 7],
[2, 3, 6, 7],
[1, 2, 5, 6],
[4, 5, 6, 7],
[1, 3, 4, 6],
[2, 3, 4, 5],
[1, 2, 4, 7]]
ddiff_prob_table_forward = [
[0, 0, 0, 0, 0, 0, 0, 0],
[1, 1, 5, 5, 3, 3, 7, 7],
[3, 7, 3, 7, 2, 6, 2, 6],
[6, 2, 2, 6, 5, 1, 1, 5],
[7, 5, 6, 4, 7, 5, 6, 4],
[4, 6, 1, 3, 6, 4, 3, 1],
[5, 3, 4, 2, 4, 2, 5, 3],
[2, 4, 7, 1, 1, 7, 4, 2]]
ddiff_prob_table_backward = [
[0, 0, 0, 0, 0, 0, 0, 0],
[1, 1, 5, 5, 3, 3, 7, 7],
[7, 3, 7, 3, 6, 2, 6, 2],
[2, 6, 6, 2, 1, 5, 5, 1],
[5, 7, 4, 6, 5, 7, 4, 6],
[6, 4, 3, 1, 4, 6, 1, 3],
[3, 5, 2, 4, 2, 4, 3, 5],
[4, 2, 1, 7, 7, 1, 2, 4]]
