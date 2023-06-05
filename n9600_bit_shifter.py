def doit(config, audio):
    audio = audio >> int(config["bit_shift"])
    return audio