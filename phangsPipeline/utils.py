"""
General purpose utilities.
"""

def select_from_list(
    master_list, 
    first=None,
    last=None,
    skip=[],
    only=[],
    loose=True,
    ):
    """
    Select only part of a list.
    """

    sorted_list = sorted(master_list,key=lambda s: s.lower())

    sub_list = []

    if first is not None:
        before_first = True
    else:
        before_first = False

    if last is not None:
        after_last = False
    else:
        after_last = False

    for element in sorted_list:

        if first is not None:
            if loose:
                if element.lower() >= first.lower():
                    before_first = False
            else:
                if element.lower() == first.lower():
                    before_first = False

        if last is not None:
            if loose:
                if element.lower() > last.lower():
                    after_last = True

        if before_first:
            continue

        if after_last:
            continue

        if skip is not None:
            if len(skip) > 0:
                match = False
                if loose:
                    for this_only in only:
                        if element.lower() == this_only.lower():
                            match = True                        
                else:
                    if element in only:
                        match = True
                if match:
                    continue

        if only is not None:
            if len(only) > 0:
                match = False
                if loose:
                    for this_only in only:
                        if element.lower() == this_only.lower():
                            match = True                        
                else:
                    if element in only:
                        match = True
                if not match:
                    continue
            
        sub_list.append(element)

        if last is not None:
            if not loose:
                if element.lower() == last.lower():
                    after_last = True

    return(sub_list)
    
