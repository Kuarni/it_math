from argparse import *
from SKUFlib.SKUF import SKUF
from SKUFlib.SVD import SVD
from SKUFlib.PowerMethodSVD import PowerMethodSVD

if __name__ == '__main__':
    parser = ArgumentParser(
        description="compress images (only BMP 24 bit) with different SVD algorithms to SKUF (Small Knowledge Unified "
                    "Format) file")

    subparsers = parser.add_subparsers(help='SKUF SKUFlib methods')


    def compress(args):
        output = args.output
        if args.output is None:
            output = args.input.rsplit('.')[0] + '.scuf'
        print("result will be saved in", output)
        algos = {
            "pySVD": SVD(compression_degree=args.compression_degree),
            "pwmSVD": PowerMethodSVD(compression_degree=args.compression_degree),
        }
        if args.algorithm not in algos:
            raise ValueError(
                f"Algorithm \"{args.algorithm}\" not supported\nList of supported algorithms: {list(algos.keys())}")
        SKUF(args.input, algos[args.algorithm]).save(output)


    p = subparsers.add_parser('compress', help="compress image to skuf file by SVD algo")
    p.add_argument('-i', "--input", help="the input file to compress", type=str, required=True)
    p.add_argument('-o', "--output", help="the output skuf file", type=str, required=False)
    p.add_argument('-a', "--algorithm",
                   help="compression SVD algorithm, there are pySVD (np linalg SVD), pwmSVD (power method SVD)",
                   type=str, required=False, default="pySVD")
    p.add_argument('-c', "--compression-degree", help="compression degree", type=float, required=False, default=2)
    p.set_defaults(func=compress)


    def decompress(args):
        output = args.output
        if args.output is None:
            output = args.input.rsplit('.')[0] + '_decompress.bmp'
        print("result will be saved in", output)
        (SKUF(args.input).decode()).save(output)


    p = subparsers.add_parser('decompress', help="decompress image from skuf file")
    p.add_argument('-i', "--input", help="the input scuf file to decompress", type=str, required=True)
    p.add_argument('-o', "--output", help="the output decompress image path", type=str, required=False)
    p.set_defaults(func=decompress)

    args = parser.parse_args()
    args.func(args)
