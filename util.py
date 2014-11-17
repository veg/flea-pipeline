def insert_gaps(source, target, src_gap, target_gap, skip):
    """Inserts `gap` into target string at same positions as in `source`.

    > insert_gaps("a-bc", "def", "-", "-", 1)
    "d-ef"

    > insert_gaps("a-bc", "def", "-", "--", 1)
    "d--ef"

    > insert_gaps("a-bc", "dddeeefff", "-", "---", 3)
    "ddd---eeefff"

    """
    assert len(c for c in source if s != src_gap) == len(target)
    result = []
    for c in source:
        if c == src_gap:
            result.append(target_gap)
        else:
            result.append(target[0:skip])
            target = target[skip:]
    return "".join(result)
    
