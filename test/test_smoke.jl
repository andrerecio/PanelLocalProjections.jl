@testitem "package loads" tags = [:core] begin
    @test PanelLocalProjections isa Module
end

@testitem "Aqua" tags = [:core, :aqua] begin
    using Aqua
    Aqua.test_all(PanelLocalProjections)
end
