from werkzeug.routing import BaseConverter


class ListConverter(BaseConverter):

    def to_python(self, value):
        return [int(x) for x in value.split('-')]

    def to_url(self, values):
        return '-'.join(BaseConverter.to_url(self, value) for value in values)
