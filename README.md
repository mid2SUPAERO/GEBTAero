# Foobar ![CI status](https://img.shields.io/badge/build-passing-brightgreen.svg)

Foobar is a Python library for dealing with word pluralization.

## Installation

### Requirements
* Linux
* Python 3.3 and up

`$ pip install foobar`

## Usage

```python
import foobar

foobar.pluralize('word') # returns 'words'
foobar.pluralize('goose') # returns 'geese'
foobar.singularize('phenomena') # returns 'phenomenon'
```

## Development
```
$ virtualenv foobar
$ . foobar/bin/activate
$ pip install -e .
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)