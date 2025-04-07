#/usr/bin/bash

CONFIG=CMakeUserPresets.json
if [ ! -f $CONFIG ]; then
echo "Generating CMakeUserPresets.json"
cat > $CONFIG <<- EOM
{
"version": 2,
"configurePresets": [
    {
    "name": "default",
    "inherits": "vcpkg",
    "environment": {
        "VCPKG_ROOT": "$VCPKG_ROOT"
    }
    }
]
}
EOM
else
echo "CMakeUserPresets.json already exists."
fi

cmake --preset=default
cmake --build build