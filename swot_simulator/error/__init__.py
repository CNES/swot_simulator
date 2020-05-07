from .. import settings


class ErrorStat:
    def __init__(self, parameters: settings.Parameters) -> None:
        self.delta_al = parameters.delta_al
        self.delta_ac = parameters.delta_ac
        self.height = parameters.height
        self.len_repeat = parameters.len_repeat
        self.nseed = parameters.nseed
