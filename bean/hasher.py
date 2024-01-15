import math


def function(hash_function_type, fragment, hash_index):
    """
    Obtain hash value from inform and hash function index by the designated type of hash function.

    :param hash_function_type: type of hash function.
    :type: string

    :param fragment: one requested information fragment.
    :type: string

    :param hash_index: index of hash function.
    :type: int

    :return: hash value.
    """
    # function directory
    # if there are some function added in this file, please add them here
    hash_func = {
        'BKDR': BKDR,
        'AP': AP,
        'DJB': DJB,
        'JS': JS,
        'SDBM': SDBM,
        'RS': RS,
        'SHA32': SHA32,
        'SHA31': SHA31,
        'nSHA31': nSHA31,
        'Murmur': Murmur,
        'Lookup3': Lookup3,
        'FNV1a': FNV_1a,
        'FNV1a_64': FNV_1a_64
    }

    if hash_function_type != 'nSHA31':
        if type(fragment[0]) is not int:
            temp = []
            for data in fragment:
                temp.append(ord(data))
            fragment = temp

    # check validity of input function type
    if hash_function_type in hash_func:
        return hash_func[hash_function_type](fragment, hash_index)

    else:
        raise Exception("Unknown hash function type %s\nValid hash function types are: %s\n" %
                        (hash_function_type, 'BKDR, AP, DJB, JS, SDBM, RS'))


def list_filling(lst):
    max_length = len(max([lst[x][y] for x in range(len(lst)) for y in range(len(lst[0]))], key=len))
    length = []
    for d1 in range(len(lst)):
        for d2 in range(len(lst[d1])):
            l = lst[d1][d2]
            length.append(len(l))
            if len(l) < max_length:
                if type(l[0]) == int:
                    lst[d1][d2] = l + [0] * (max_length - len(l))
                else:
                    lst[d1][d2] = l + 'A' * (max_length - len(l))

    return lst, length


# noinspection PyPep8Naming
def BKDR(fragment, hash_index):
    """
    Obtain hash value by BKDR Hash.

    :param fragment: one requested inform.
    :type: string

    :param hash_index: index of hash function.
    :type: int

    :return: hash value.
    """
    seed = 131
    value = hash_index

    for data in fragment:
        value = value * seed + data

    return value & 0x7FFFFFF


# noinspection PyPep8Naming
def AP(fragment, hash_index):
    """
    Obtain hash value by AP Hash.
    """
    value = hash_index

    for index, data in enumerate(fragment):
        if index % 2 == 0:
            value ^= (value << 7) + data + (value >> 3)
        else:
            value ^= (~(value << 11)) + data + (value >> 5)

    return value & 0x7FFFFFF


# noinspection PyPep8Naming
def DJB(fragment, hash_index):
    """
    Obtain hash value by DJB Hash.
    """
    value = 5381 + hash_index

    for data in fragment:
        value += (value << 5) + data

    return value & 0x7FFFFFF


# noinspection PyPep8Naming
def JS(fragment, hash_index):
    """
    Obtain hash value by JS Hash.
    """

    value = 1315423911 + hash_index

    for data in fragment:
        value ^= (value << 5) + data + (value >> 2)

    return value & 0x7FFFFFF


# noinspection PyPep8Naming
def RS(fragment, hash_index):
    """
    Obtain hash value by RS Hash.
    """

    base_1 = 63689
    base_2 = 378551
    value = hash_index

    for data in fragment:
        value = value * base_1 + data
        base_1 *= base_2

    return value & 0x7FFFFFF


# noinspection PyPep8Naming
def SDBM(fragment, hash_index):
    """
    Obtain hash value by SDBM Hash.
    """
    value = hash_index

    for data in fragment:
        value = data + (value << 6) + (value << 16) - value

    return value & 0x7FFFFFF


# noinspection PyPep8Naming
def SHA31(fragment, hash_index):
    """
    Obtain hash value by SHA31 Hash.
    """
    hash_value = 12589648753214523651 + hash_index

    for data in fragment:
        hash_value *= 1254145215421
        hash_value ^= data
        hash_value %= 12589648753214523651

    return hash_value


def DNA2Number(seq):
    """
    Convert a DNA fragment(less than 32nt) to integer.
    """
    rz = 0
    base2num = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    for i in range(-1, len(seq)*-1-1, -1):
        rz += base2num[seq[i]] * math.pow(4, -1-i)

    return int(rz)


# new version of SHA31 called nSHA31
def nSHA31(fragment, hash_index):
    """
    Obtain hash value by nSHA31 Hash.
    """

    base_value = 12589648753214523651
    times_value = 1099511628211
    device_value = 14695981039346656037
    hash_value = base_value + hash_index

    for i in range(0, len(fragment), 32):
        hash_value *= times_value
        hash_value ^= DNA2Number(fragment[i:i+32])
        hash_value %= device_value

    return hash_value


def SHA32(fragment, hash_index):
    """
    Obtain hash value by SHA32 Hash.
    """
    hash_value = 1099511628211 + hash_index

    for data in fragment:
        hash_value *= 1099511628211
        hash_value ^= data
        hash_value %= 14695981039346656037

    return hash_value


def Murmur(fragment, hash_index):
    """
    Obtain hash value by Murmur3_32bit Hash.
    """
    c1 = 0xcc9e2d51
    c2 = 0x1b873593
    r1 = 15
    r2 = 13

    fragment += [0 for _ in range(4 - len(fragment) % 4)]

    for i in range(0, len(fragment), 4):
        # little endian load order
        k = (fragment[i + 3] << 24) ^ (fragment[i + 2] << 16) ^ (fragment[i + 1] << 8) ^ fragment[i]
        k *= c1
        k %= 4294967296
        k = (k << r1) | (k >> 32 - r1)  # ROTL32(k1,r1)
        k *= c2
        k = k % 4294967296
        hash_index ^= k
        hash_index = (hash_index << r2) | (hash_index >> 32 - r2)  # ROTL32(h1,r2)
        hash_index = hash_index * 5 + 0xe6546b64
        k %= 4294967296

    k = 0

    if (len(fragment) & 0x03) == 3:
        k = fragment[int(len(fragment) / 4 + 2)] << 16
    if (len(fragment) & 0x03) >= 2:
        k ^= fragment[int(len(fragment) / 4 + 1)] << 8
    if (len(fragment) & 0x03) >= 1:
        k ^= fragment[int(len(fragment) / 4)]
        k *= c1
        k = k % 4294967296
        k = (k << r1) | (k >> 32 - r1)  # ROTL32(k1,15);
        k *= c2
        k = k % 4294967296
        hash_index ^= k

    hash_index ^= len(fragment)
    hash_index ^= hash_index >> 16
    hash_index *= 0x85ebca6b
    hash_index = hash_index % 4294967296
    hash_index ^= hash_index >> 13
    hash_index *= 0xc2b2ae35
    hash_index = hash_index % 4294967296
    hash_index ^= hash_index >> 16

    return hash_index % 4294967296


def Lookup3(fragment, hash_index):
    """
    Obtain hash value by Lookup3 Hash.
    """

    hash_value = hash_index

    for data in fragment:
        hash_value += data
        hash_value += (hash_value << 10)
        hash_value ^= (hash_value >> 6)

    hash_value += (hash_value << 3)
    hash_value ^= (hash_value >> 11)
    hash_value += (hash_value << 15)

    return hash_value & 0x7FFFFFF


def FNV_1a(fragment, hash_index):
    """
    Obtain hash value by FNV1a_32bit Hash.
    """
    hash_value = 2166136261 + hash_index

    for data in fragment:
        hash_value ^= data
        hash_value *= 16777619

    return hash_value & 0x7FFFFFF


def FNV_1a_64(fragment, hash_index):
    """
    Obtain hash value by FNV1a_64bit Hash.
    """
    hash_value = 14695981039346656037 + hash_index

    for data in fragment:
        hash_value ^= data
        hash_value *= 1099511628211

    return hash_value & 0x7FFFFFF
