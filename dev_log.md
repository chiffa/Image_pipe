# This is a general library for image processing

## Logic in the 

Logic layers:
  1. Traversal
  2. Matching of the images to different channels
  3. Segmentation of
      a. Mitochondria
      b. Cells
  4. Retrieval of descriptors each surface/volume
  5. Classification of a surface/volume
      a. with respect to the signal intensity => Cell death
      b. with respect to the geometry
          i. Area, skeleton, =>  fragmented mitochondria or not
  6. Re-hashing of the file names in order to extract time-stamp
  7. Perform time-line tracing and statistics calculation

  DEBUG: output the graphical representation of what has been done by different cell


## Desired code:

```python
# Attach the debuger
set.debugger(dump_location = "./verification/", layers={"layer1":12, })

# point to the source folder 
set.source_folder("xxxx")

# map a pattern in the name to a color channel
input.images = set.color_pattern({"_1465":"green", "_2660":"red"})

intermediate.cells = pipeline.step.add(segment.cells(input.images["green"], binary=True), debugger=layer1)
intermediate.cells = pipeline.step.add(remove outliers(intermediate.cells, "green"), debugger=layer2)
intermediate.cells = pipeline.step.add(remove outliers(intermediate.cells, "red"), debugger=layer2)
intermediate.granules = pipeline.step.add(segment.granules(intermediate.cells, "green"), debugger=layer3)
intermediate.mitochondriae = pipeline.step.add(segment.threshold(intermediate.cells, "red"), debugger=layer3)

colocalization = pipeline.step.add(colocalization.first_in_second(granules, mitochondriae), debugger=layer3)
anti_colocalization = pipeline.step.add(colocalization.first_in_second(granules, not(mitochondriae)), debugger=layer3)

pipeline.save.pattern(input.image.shorcode, colocalization, anti_colocalization)

save.location("final_table.csv")

pipeline.run()
```

## Remarks:

In fact, we don't really need an explicit pipeline class to assem ble the pipeline.
We use intermediate variables within the user class in order `to stitch the pipeline 
and store the intermediate processes. We can as well plug the debug process to 
get the function debug output as it goes.

Unless we are using `yield` statements, in which case the pipeline would be performing
processing ticks.

Idea - separate the namespaces that are logical in the context of the analysis

In order to avoid the generator once-consumed property, use multiple pass-through 
objects in dicts. That will also make it for an easy debugging in the end