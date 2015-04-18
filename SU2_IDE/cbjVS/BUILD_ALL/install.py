import os 
import optparse
import shutil

def main():
    parser = optparse.OptionParser()
    parser.add_option("-s",dest="source_dir")
    parser.add_option("-d",dest="dest_dir")
    options, args = parser.parse_args()
    install_binary(options.source_dir,options.dest_dir)

def install_binary(dir_from, dir_to):
    files = ['CBJ_DEF','SU2_CFD','SU2_DEF','SU2_DOT','SU2_GEO','SU2_MSH','SU2_SOL']
    for fn in files:
        source_file = os.path.join(dir_from,fn+'.exe')
        dest_file = os.path.join(dir_to,fn+'.exe')
        if not os.path.exists(dir_to):
            os.makedirs(dir_to)
        if os.path.exists(source_file):
            if os.path.exists(dest_file):
                os.remove(dest_file)
            shutil.copy(source_file,dest_file)
        else:
            print 'File:{} not found!'.format(fn)

if __name__=="__main__":
    main()

