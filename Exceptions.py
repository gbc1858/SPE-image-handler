class CurrentValueError(BaseException):
    def __init__(self, message):
        super().__init__(self, message)


class DenoiseBGError(BaseException):
    def __init__(self, message):
        super().__init__(self, message)


class DenoiseCroppedError(BaseException):
    def __init__(self, message):
        super().__init__(self, message)


class DenoiseFrameAfterContourError(BaseException):
    def __init__(self, message):
        super().__init__(self, message)


class RmsMethodError(BaseException):
    def __init__(self, message):
        super().__init__(self, message)