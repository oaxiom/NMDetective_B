
import random
import matplotlib.pyplot as plot

def load_tid_table(filename):
    oh = open(filename, 'rt')

    res = {}

    for line in oh:
        line = line.strip().split('\t')
        res[line[0]] = float(line[1])
    oh.close()

    return res



obs = load_tid_table('observed.tidtable.tsv')
exp = load_tid_table('expected.Lindeboom.tidtable.tsv')

# match up keys:
matching_keys = list(set(obs.keys() & exp.keys()))

print(f'Observed = {len(obs)}')
print(f'Expected = {len(exp)}')
print(f'Matching = {len(matching_keys)}')

m = {}
for k in matching_keys:
    m[k] = (obs[k], exp[k])

fig = plot.figure()
ax = fig.add_subplot(111)

x = [m[k][0] for k in sorted(m.keys())]
y = [m[k][1] for k in sorted(m.keys())]

s = [random.randint(10,2000) for i in range(len(x))]

ax.scatter(x, y, s=s, alpha=0.01, ec='none')
ax.set_xlabel('Observed')
ax.set_ylabel('Expected')

fig.savefig('scatter.pdf')

