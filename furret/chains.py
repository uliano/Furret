class SimpleRange:  # contiuguous range of residues (extremes included)

    @staticmethod
    def is_valid_range(the_range):
        if len(the_range) != 2:
            return False
        begin, end = the_range
        if 0 < begin <= end:
            return True

    def __repr__(self):
        return f"SimpleRange({self.begin}, {self.end})"

    def __init__(self, the_range):
        if SimpleRange.is_valid_range(the_range):
            self.begin = int(the_range[0])
            self.end = int(the_range[1])
        else:
            raise ValueError(f"Invalid range(s) {the_range} in SimpleRange.__init__()")

    def __gt__(self, other):  # self follows other
        assert isinstance(other, SimpleRange)
        if self.begin == other.begin:
            return self.end > other.end
        else:
            return self.begin > other.begin

    def __ge__(self, other):  # self follows other
        assert isinstance(other, SimpleRange)
        if self.begin == other.begin:
            return self.end >= other.end
        else:
            return self.begin >= other.begin

    def __lt__(self, other):  # self precedes other
        assert isinstance(other, SimpleRange)
        if self.begin == other.begin:
            return self.end < other.end
        else:
            return self.begin < other.begin

    def __le__(self, other):  # self precedes other
        assert isinstance(other, SimpleRange)
        if self.begin == other.begin:
            return self.end <= other.end
        else:
            return self.begin <= other.begin

    def __eq__(self, other):
        assert isinstance(other, SimpleRange)
        return self.begin == other.begin and self.end == other.end

    def __iter__(self):
        for i in range(self.begin, self.end + 1):
            yield i

    def __contains__(self, other):
        assert isinstance(other, SimpleRange)
        return self.begin <= other.begin and self.end >= other.end

    def __len__(self):
        return self.end - self.begin + 1

    def union(self, other):
        assert isinstance(other, SimpleRange)
        # check for intersection points or contiguity
        if (self <= other and other.begin - self.end <= 1) or (self > other and self.begin - other.end <= 1):
            begin = min(self.begin, other.begin)
            end = max(self.end, other.end)
            the_range = (begin, end)
            try:
                return [SimpleRange(the_range)]
            except ValueError:
                return []
        else:
            return sorted([self, other])


class Chain:  # a named list of SimpleRanges in ascending order
    def __init__(self, name, ranges):
        self.name = name
        self.ranges = ranges if type(ranges) is list else [ranges]
        self._reduce_ranges()

    def _reduce_ranges(self):
        if self.ranges:
            self.ranges.sort()
            reduced_ranges = [self.ranges[0]]
            for r in self.ranges[1:]:
                reduced_ranges += r.union(reduced_ranges.pop())
            self.ranges = reduced_ranges

    def __eq__(self, other):
        return other.name == self.name and other.ranges == self.ranges

    def __repr__(self):
        return f'Chain {self.name}({self.ranges})'

    def __len__(self):
        return sum(len(r) for r in self.ranges)

    def add_ranges(self, ranges):
        ranges = ranges if type(ranges) is list else [ranges]
        self.ranges += ranges
        self._reduce_ranges()

    def union(self, other):
        assert isinstance(other, Chain)
        if self.name == other.name:
            me = Chain(self.name, self.ranges + other.ranges)
            return [me]
        else:
            return [self, other]


class ChainGroups:
    """
    Class for manipulating groups of chains.

    init:list of strings like ["A/B=96-516", "C/D/E=9-94" ...]

    supported methods:
    chain_names(self) -> generator returning chain names
    get_chain(self, name) -> return Chain object by its name
    add_chain(self, other) -> add a single Chain object
    add_chains_from_string(self, string) -> add one or multiple chains from a single string like 'A/B/B=1-100'
    __iter__(self) --> yield single Chain objects
    """

    def __init__(self, chain_group_list):  # chain group list ["A/B=96-516","C/D/E=9-94"]
        self.chains = []
        for element in chain_group_list:
            self.add_chains_from_string(element)

    def chain_names(self):
        for c in self.chains:
            yield c.name

    def get_chain(self, name):
        for c in self.chains:
            if c.name == name:
                return c

    def add_chain(self, other):
        assert isinstance(other, Chain)
        if other.name in self.chain_names():
            self.get_chain(other.name).add_ranges(other.ranges)
        else:
            self.chains.append(other)

    def add_chains_from_string(self, string):
        key, values = string.split('=')
        chains = key.split('/')
        residue_range = SimpleRange(tuple(map(lambda x: int(x), values.split('-'))))
        for chain_name in chains:
            if chain_name in self.chain_names():
                self.get_chain(chain_name).add_ranges(residue_range)
            else:
                self.chains.append(Chain(chain_name, residue_range))

    def __repr__(self):
        return f"GhainGroups({self.chains})"

    def __iter__(self):
        for chain in self.chains:
            yield chain

    def merged(self):
        if not self.chains:
            return None
        merged_ranges = []
        chain_names = []
        for chain in self.chains:
            merged_ranges += chain.ranges
            chain_names.append(chain.name)
        return Chain(''.join(chain_names), merged_ranges)
