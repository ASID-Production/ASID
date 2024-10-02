from rest_framework import mixins, viewsets


class StructureModelViewSet(
    # mixins.CreateModelMixin,
    mixins.RetrieveModelMixin,
    # mixins.UpdateModelMixin,
    mixins.DestroyModelMixin,
    mixins.ListModelMixin,
    viewsets.GenericViewSet
):
    '''
    A viewset that provides default `retrieve()`, `destroy()` and `list()` actions.
    '''

    pass
