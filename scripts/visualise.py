import stratiMesh
import argparse

def main():

    parser = argparse.ArgumentParser(
        prog="visualise",
        description="""Convert outputs from Carst to 3D visualisation for paraview"""
    )
    #Parse command-line argumens.
    parser.add_argument('output_directory', 
                        type=str, 
                        help='Name of the output directory used in Carst'
                        )
    # TODO: we could make this either the options file or a directory and check which it is. If options file, parse and pull
    # out the output directory. However, this might be tricky given we would need a full path. This would make the
    # next few options redudant...
    parser.add_argument('start_time', 
                           type=float, 
                           help='Start time of the simulation'
                        )
    parser.add_argument('end_time', 
                        type=float, 
                        help='End time of the simulation.'
                        )
    parser.add_argument('dispTime',
                        type=float,
                        help="How often if your Carst model output? Defaults to 500 years",
                        )
    parser.add_argument('-v', 
                        dest='verbose', 
                        default=False, 
                        action='store_true',
                        help='verbose output'
                        )
    parser.add_argument('--name',
                        '-n',
                        dest='xdmfName', 
                        default='stratal_series',
                        help='Output xdmf name. Default is stratal_series'
                        )
    
    args = parser.parse_args()
    folder = args.output_directory
    xdmfName = args.xdmfName
    verbose = args.verbose
    dispTime = args.dispTime
    startTime = args.start_time
    endTime = args.end_time

    mesh = stratiMesh.stratiMesh(folder = folder, xdmfName = xdmfName, dispTime = dispTime, verbose=verbose)
    mesh.outputSteps(startTime = startTime, endTime = endTime)


if __name__=='__main__':
    main()
