import base64
import os
from cryptography.fernet import Fernet
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC


def hash_password(password):
    """
    This is not secure against dictionary attacks as I am not sharing the salt.
    """
    salt=b'\xccn\x8c\xcaj-\xb8\x8fL&\x15\x04\x8bo\xc6/'
    password = password.encode()
    kdf = PBKDF2HMAC(
        algorithm=hashes.SHA256(),
        length=32,
        salt=salt,
        iterations=390000,
    )
    return base64.urlsafe_b64encode(kdf.derive(password))

def encrypt(message: str, password: str):
    key = hash_password(password)
    f = Fernet(key)
    return f.encrypt(message.encode())

def decrypt(token: bytes, password: str):
    key = hash_password(password)
    return Fernet(key).decrypt(token).decode()
