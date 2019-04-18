

T = 1.0
nT = 10000.0


def test(t, h1, h2):
    d = t/pow(h1, 2) + t/pow(h2, 2)
    return d <= 0.5


for x in range(10, 1000, 10):
    if test(T/nT, 1./x, 1./x):
        print("!!!", x)
