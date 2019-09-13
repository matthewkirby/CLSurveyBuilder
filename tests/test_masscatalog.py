import clsurveybuilder as cls

if __name__ == "__main__":
    builder = cls.SurveyBuilder(None)
    mock = builder.build_survey(None, 1, save_m200=True)
    print(mock)
