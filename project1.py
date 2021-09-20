import math
import random

def generate_pseudo_primes():
    """
    The generate_pseudo_primes function uses Fermat's little theorem to find
    two large pseudo prime numbers

    Returns
    -------
    primes : list
        A list containing two prime numbers.

    """
    primes = [] # Empty list to hold our two primes
    while len(primes) != 2:
        candidate = random.randint(10000, 50000) # Generate number to test
        composite = False

        # Run Fermat's test 20 times
        # If all 20 tests pass, then the probability the candidate is not prime is less than 1/2^20,
        # and we add it to our list
        for i in range(20):
            # Randomly generate an a in (2, n - 1), then check if a^(n-1) mod n = 1, where n = candidate
            if pow(random.randint(2, candidate - 1), candidate - 1, candidate) != 1:
                composite = True
        if not composite:
	        primes.append(candidate)
    return primes

def extended_gcd(a, b):
    """ 
    The extended_gcd function implements the
    extension of Euclid's GCD algorithm to find integers x and y
    such that ax + by = gcd(a, b) 
    input: Positive Integers a, b such that a = e and b = phi
    output: Integers x, y, d such that d = gcd(a, b) and ax + by = d
    @author: Dr. Hu
    """
    if b == 0:
        return (1, 0, a)
    (x, y, d) = extended_gcd(b, a % b)
    return y, x - a // b * y, d

def generate_keys():
    """
    The generate_keys function generates the public and private keys

    Returns
    -------
    list
        A list of length two that contains the public key (n, e) and the 
        private key d.

    """
    # step 1
    # find two large prime numbers using Fermat's little theorem
    p, q = generate_pseudo_primes()
    
    # Calculate n and phi
    n = p * q
    phi = (p - 1) * (q - 1)

    # step 2
    # select a public key e
    # e should be relatively prime to phi
    # that is gcd(e, phi) = 1
    # Randomly pick e in (2, phi - 1)
    found_e = False
    while not found_e:
        x = random.randint(2, phi - 1)
        if math.gcd(x, phi) == 1:
            e = x
            found_e = True    

    # step 3
    # determine private key d
    # d is multiplicative inverse of e in Zphi
    # that is (e * d) mod phi = 1
    # to find d use extended Euclid's algorithm
    # input a = e and b = phi
    # We apply the extended Euclidean algorithm to find the x and y
    # in the equation 1 = e*x + phi*y.
    # The x is the private key d because phi*y mod phi = 0.
    d = extended_gcd(e, phi)[0] % phi
    
    # Return the public key(n, e) and private key d
    return [(n, e), d]

def encrypt(plaintext, public_key):
    """
    The encrypt function encrypts the plaintext using the public key.

    Parameters
    ----------
    plaintext : string
        The message to be encrypted.
    public_key : tuple
        The public key (n, e).

    Returns
    -------
    cyphertext : list
        A list of integers representing each encrypted character of the plaintext.

    """
    # Convert each character of plaintext to an integer 
    # Then apply encryption formula C = M^e mod n using fast modulo exponentiation algorithm
    # Store each encrypted character in a list
    # Return the list
    cyphertext = [pow(ord(ch), public_key[1], public_key[0]) for ch in plaintext]
    return cyphertext

def decrypt(cyphertext, d, n):
    """
    The decrypt function decrypts the cyphertext using the private key d and the modulus n.

    Parameters
    ----------
    cyphertext : list
        A list of integers representing each encrypted character of the plaintext.
    d : integer
        The private key.
    n : integer
        The modulus.

    Returns
    -------
    plaintext : string
        The original plaintext.

    """
    # For each integer in list of cyphertext, apply the decryption formula M = C^d mod n using fast modulo exponentiation algorithm
    # Then convert resulting integer to ascii
    plaintext = ''
    for ch in cyphertext:
        plaintext += chr(pow(ch, d, n))
    return plaintext

def sign(M, d, n):
    """
    The sign function takes the message, private key, and the modulus and generates
    a signature. The signature = m^d mod n.

    Parameters
    ----------
    M : string
        The message.
    d : int
        Private key.
    n : int
        The modulus n (n=p*q).

    Returns
    -------
    signature : list
        The signature is a list of integers.

    """
    # Signature = M^d mod n
    signature = [pow(ord(ch), d, n) for ch in M]
    return signature


def verify(M, e, n, s):
    """
    The verify function takes the sent message, the signature, the public key e, 
    and the modulus n and produces a message. That message is compared against 
    the original message for verification. It returns a boolean result of that 
    comparison.

    Parameters
    ----------
    M : string
        The message sent with the signature.
    e : int
        The public key.
    n : int
        The modulus.
    s : list
        The signature.

    Returns
    -------
    bool
        A boolean that represents whether the message sent with the signature
        matches the original message.

    """
    # Message_prime = Signature^e mod n
    m_prime = '' # string to hold message built from undoing the signature
    for ch in s:
        m_prime += chr(pow(ch, e, n))
    return M == m_prime

def get_input():
    restart = 'Y'
    while restart == "Y" or restart == "YES":
        m = ''
        public_key, private_key = generate_keys()
        print()
        # Prompt user for message
        # If message is null or empty, Re-prompt user
        # Else continue and encrypt message
        while (len(m.strip()) == 0):
            m = input('Enter your message to be encrypted: ') #Alice's message to Bob
        cyphertext = encrypt(m, public_key)
        print(f'\nPublic Key: {public_key}')
        print(f'Private Key: {private_key}')
        print(f'Encrypted message: {cyphertext}') #What Trudy should see
        message = decrypt(cyphertext, private_key, public_key[0])
        print(f'\nMessage decrypted: {message}') #What Bob receieves from Alice
        signature = sign(m, private_key, public_key[0])
        print(f'\nSignature made using private key: {signature}')
        verification = verify(m, public_key[1], public_key[0], signature)
        print(f'\nDoes verifing the signature match the orginal message: {verification}')
        # Prompt user to another message
        # If 'Yes' prompt for message
        # Else end program
        restart = (input('\nSend another message (Y/N)?: ')).upper()
        if restart == "Y" or restart == "YES":
            print('\n------------------------------------------')
        else:
            break

get_input()
    


