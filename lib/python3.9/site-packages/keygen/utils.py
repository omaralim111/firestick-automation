import string


def build_charset(signature: str) -> str:
    charset = ''
    if 'u' in signature:
        charset += string.ascii_uppercase
    if 'l' in signature:
        charset += string.ascii_lowercase
    if 'd' in signature:
        charset += string.digits
    if 'p' in signature:
        charset += string.punctuation
    if 'h' in signature:
        charset += string.hexdigits
    if 'o' in signature:
        charset += string.octdigits
    if 'b' in signature:
        charset += '01'

    return ''.join(sorted(set(charset)))
