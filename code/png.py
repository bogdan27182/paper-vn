import zlib, struct
PAL = {0: (255, 255, 255), 1: (220, 30, 30), 2: (40, 40, 220)}
def write_png(fn, rows, pal=PAL, scale=1):
    h = len(rows); w = len(rows[0])
    raw = bytearray()
    for r in rows:
        line = bytearray([0])
        for v in r:
            line += bytes(pal[v]) * scale
        for _ in range(scale):
            raw += line
    def chunk(t, d):
        c = struct.pack('>I', len(d)) + t + d
        return c + struct.pack('>I', zlib.crc32(t + d) & 0xffffffff)
    png = b'\x89PNG\r\n\x1a\n' + chunk(b'IHDR', struct.pack('>IIBBBBB', w*scale, h*scale, 8, 2, 0, 0, 0))
    png += chunk(b'IDAT', zlib.compress(bytes(raw), 6)) + chunk(b'IEND', b'')
    open(fn, 'wb').write(png)
