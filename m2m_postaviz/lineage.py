class Lineage():
    def __init__(self):
        self.level_dict=dict()

    def get_keys(self, data:dict, level:int):
        keys = data.keys()
        keys_with_child = []
        for k in keys:
            if bool(data[k]):
                keys_with_child.append(k)
        if level in self.level_dict:
            self.level_dict[level].extend(keys_with_child)
        else:
            self.level_dict[level] = keys_with_child

    def construct_dict(self, data:dict, level:int):
        for k in data:
            self.get_keys(data[k], level)
            self.construct_dict(data[k], level+1)