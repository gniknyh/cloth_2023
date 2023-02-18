# import pyweaving.wif

from pyweaving import Draft, instructions
from pyweaving.wif import WIFReader, WIFWriter
from pyweaving.render import ImageRenderer


pixel_scale = 4

def load_draft(infile):
    if infile.endswith('.wif'):
        return WIFReader(infile).read()
    elif infile.endswith('.json'):
        with open(infile) as f:
            return Draft.from_json(f.read())
    else:
        raise ValueError(
            "filename %r unrecognized: .wif and .json are supported" %
            infile)


# def render(opts):
#     draft = load_draft(opts.infile)
#     if opts.outfile:
#         if opts.outfile.endswith('.svg'):
#             SVGRenderer(draft).save(opts.outfile)
#         else:
#             ImageRenderer(draft).save(opts.outfile)
#     else:
#         ImageRenderer(draft).show()


def render(file):
    draft = load_draft(file)
    ir = ImageRenderer(draft)
    ir.set_pixel_resolution(pixel_scale)
    # ir.show()
    ir.save('out.png')


def convert(opts):
    draft = load_draft(opts.infile)
    if opts.outfile.endswith('.wif'):
        WIFWriter(draft).write(opts.outfile)
    elif opts.outfile.endswith('.json'):
        with open(opts.outfile, 'w') as f:
            f.write(draft.to_json())


def thread(opts):
    draft = load_draft(opts.infile)
    instructions.threading(draft, opts.repeats)


def weave(opts):
    draft = load_draft(opts.infile)
    assert opts.liftplan, "only liftplan supported for now"
    save_filename = '.' + opts.infile + '.save'
    print("SAVE FILENAME is %r" % save_filename)
    instructions.weaving(draft,
                         repeats=opts.repeats,
                         start_repeat=opts.start_repeat,
                         start_pick=opts.start_pick,
                         save_filename=save_filename)


def tieup(opts):
    draft = load_draft(opts.infile)
    instructions.tieup(draft)


def stats(opts):
    draft = load_draft(opts.infile)
    warp_longest, weft_longest = draft.compute_longest_floats()
    print("Title:", draft.title)
    print("Author:", draft.author)
    print("Address:", draft.address)
    print("Email:", draft.email)
    print("Telephone:", draft.telephone)
    print("Fax:", draft.fax)
    print("Notes:", draft.notes)
    print("Date:", draft.date)
    print("***")
    print("Warp Threads:", len(draft.warp))
    print("Weft Threads:", len(draft.weft))
    print("Shafts:", len(draft.shafts))
    print("Treadles:", len(draft.treadles))
    print("Longest Float (Warp):", warp_longest)
    print("Longest Float (Weft):", weft_longest)

render("05sH010_linen.wif")