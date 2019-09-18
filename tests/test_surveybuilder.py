# import clsurveybuilder as cls
import numpy as np
import sys
sys.path.append("../clsurveybuilder")
import surveybuilder as cls




def main():
    str2sdeg = 180.**2/np.pi**2

    builder = cls.SurveyBuilder('sample_config.ini')
    print('cosmo:', builder.cosmo)
    print('zmin:', builder.zmin)
    print('zmax:', builder.zmax)
    print('solidA:', builder.solid_angle*str2sdeg)
    print('name:', builder.surveyname)
    print('path:', builder.outputpath)

    builder.build_survey(None, 1, save_m200=True)
    return







if __name__ == "__main__":
    main()
