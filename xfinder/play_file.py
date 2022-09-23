def test_sets():

    x = ("a", 1, 2), ("b", 3, 4)
    sublists_queue = set(x)
    sublists_in_cluster = set().union(sublists_queue)

    print(sublists_queue)
    print(sublists_in_cluster)

x = 5
z = [{2},{x, 3, 4}, {0, 1}]

z = sorted(z, key=len, reverse=True)


print(z)