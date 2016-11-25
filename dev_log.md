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

## Formalism:
Can we replace the in_channel, out_channel, log_channel with actual bindings to the variables?
Right now, we are using an explicit dict for scope passing - could we rather use the function names?

What I am trying is basically to make flow control explicit by using an assembly-line 
model. The sort-coming is that I am using words, so I am using a lot Python interpreter errors
and IDEs power.

A way around it is to hack the pipeline in the same manner as the guy who wrote hy (py-lisp)

## Sub-segmentation
As of now, the mechanics of our method are operation only on the full image frames.
If we want to segment the image for the analysis any further, for instance in order to
perform cell-specific processing, would need to:
    - have a sencond generator wrapping / unwrapping routine
    - second-level dictionary store in the dict. Basically a 
    
The second generator is basically a wrapper around the generator wrapper to execute 
the same instruction over a generator of generators, where the inner generator is processed
by the generator wrapper logic, whereas the outer is processed by the second layer of wrapping

=> That's basically a double generator_wrapper, but that cannot support the expected dims because
a dict is being passed.
=> We actually just need the outer wrapper, the inner can be just the generator wrapper

## Registering Run:
As of now, it is pretty trivial to get the current revision from git and the source/dump locations
of all the elements. to get the replicability within the pipeline.

## Possible usages:
- Assembly of generator-based sub-pipes
- Assembly of input-output chains and then wrapping them into pipe

## Organization:
- Core functions => Non-wrapped, testable
- Pipeline logic => Generator-wrapped, ready for assembly; wrappers
- Pre-assembled filters
- It is all to be imported into the actual field and 

## From the usage
- we definitely need a pipeline assembly - it is too easy to forget the pipe redirections between the generators
- it is a bit frustrating to be unable to add elementary modifications to channels when they are injected
- protection against dims mismatch is good and saves some time
- the fact that stack trace returns nothing informative is definitely a minus

## Reformulation
A pretty clear way of dealing with it is to re-write to get rid of the wrappers
and reduce it all to a main for loop with embedded for sub-loops.
Nitty-gritty details:
    - splitter needs to be a generator, returning values we want to use in the end
    - secondary namespace
    - point/tile/summarize

## Audit:
- Word audit trail explaining what function does?
- Image binding into the final rendering frame?
