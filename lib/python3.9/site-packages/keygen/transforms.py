import hashlib
import uuid
from enum import Enum


class Algorythm(Enum):
    MD5 = 'md5'
    SHA1 = 'sha1'
    SHA224 = 'sha224'
    SHA256 = 'sha256'
    SHA384 = 'sha384'
    SHA512 = 'sha512'
    UUID = 'uuid'

    def __str__(self) -> str:
        return self.value

    def apply(self, value: str) -> str:
        if self in (self.MD5, self.SHA1, self.SHA224, self.SHA256, self.SHA384, self.SHA512):
            return getattr(hashlib, self.value)(value.encode()).hexdigest()
        if self is self.UUID:
            return str(uuid.uuid3(uuid.NAMESPACE_OID, value))
