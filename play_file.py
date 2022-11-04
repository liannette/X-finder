def test_sets():

    x = ("a", 1, 2), ("b", 3, 4)
    sublists_queue = set(x)
    sublists_in_cluster = set().union(sublists_queue)

    print(sublists_queue)
    print(sublists_in_cluster)

def no_return():
    if False:
        return "Never return"
    
tpl = ("ref", 123, 456)
all_sublists = [("ref", 123, 456), ("ref", 124, 456)]
dct = dict.fromkeys(all_sublists, None)
print(dct)